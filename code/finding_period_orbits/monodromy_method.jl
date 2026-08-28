import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2,
      NonlinearSolve, ADTypes, DataFrames, Dates, ProgressMeter, Printf, ColorSchemes
include("../styles/makie_theme.jl")
set_theme!(QUARTO_THEME)
include("../models/henon_heiles.jl")
using .HenonHeiles

BLAS.set_num_threads(1)          # your linear algebra is 2x2; BLAS threads only compete

const CONFIG_DIR = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
const DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
const FIG_DIR    = joinpath(@__DIR__, "../../figures/henon-heiles/periodic-orbits/")
const SAVE_DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/periodic-orbits/")

data_file        = joinpath(DATA_DIR, "E0.1127-T10000.0-py0.0-n256.jld2")

configurations = JSON3.read(read(CONFIG_DIR, String))
cfg            = configurations.explore
param          = [Float64(cfg.a.value), Float64(cfg.m.value), Float64(cfg.w.value)]
dt             = Float64(cfg.dt.value)
x0             = Float64(cfg.x0.value)

const RGRID  = range(-1.0, 1.0, length = 240)
const EPOT   = [HenonHeiles.potential(x, y, param) for x in RGRID, y in RGRID]
const LEVELS = collect(logrange(5.0*0.009, 6.9*0.089, 7))

struct SectionParams{IF, ID, P}
    integ_fast  :: IF
    yf          :: Vector{Float64}
    pyf         :: Vector{Float64}
    tsf         :: Vector{Float64}
    nmax_fast   :: Base.RefValue{Int}
 
    integ_dense :: ID
    yd          :: Vector{Float64}
    pyd         :: Vector{Float64}
    tsd         :: Vector{Float64}
    nmax_dense  :: Base.RefValue{Int}
 
    E           :: Float64
    p           :: P
    tmax        :: Float64
end



const Row = @NamedTuple begin
    E::Float64; n::Int
    seed_y::Float64; seed_py::Float64
    y::Float64; py::Float64
    prime::Int; T::Float64
    trace::Float64; detDT::Float64
    resnorm::Float64; iters::Int
    sec_y::Vector{Float64}; sec_py::Vector{Float64}
    history_y::Vector{Float64}; history_py::Vector{Float64}
    traj_x::Vector{Float64}; traj_y::Vector{Float64}
    class::String; index::Int
    id::Int; origin::String
end
 
""" 
    E::Float64; n::Int
    seed_y::Float64; seed_py::Float64
    y::Float64; py::Float64
    prime::Int; T::Float64
    trace::Float64; detDT::Float64
    resnorm::Float64; iters::Int
    sec_y::Vector{Float64}; sec_py::Vector{Float64}
    history_y::Vector{Float64}; history_py::Vector{Float64}
    traj_x::Vector{Float64}; traj_y::Vector{Float64}
    class::String; index::Int
    id::Int; origin::String
"""
orbit_table() = DataFrame(Row[])


# --- background section ---
# data = load(data_file, "results")
# y_all, py_all = Float32[], Float32[]
# for d in data
#     append!(y_all,  d.sec_y)
#     append!(py_all, d.sec_py)
# end

const EPS_OFF = 1e-9

"returns -> 1 := elliptic, 2 := hyperbollic (stable), 3 := hyperbollic (unstable), 4 := parabolic,"
function kind_index(τ; ε = 1e-6)
    abs(abs(τ) - 2) < ε && return 4   # parabolic
    abs(τ) < 2 && return 1            # elliptic
    τ > 0 ? 2 : 3                     # hyperbollic (unstable, stable)
end

const KIND_LS    = [:solid, :dash, :dashdot, :dot]
const KIND_MS    = [:circle, :xcross, :diamond, :utriangle]
const KIND_LABEL = ["elliptic", "hyperbolic", "inverse hyperbolic", "parabolic"]


const COLOR_SCHEME = [C_ORANGE, C_TEAL, C_CREAM, C_GOLD, C_PURPLE, C_GREEN]
pick_color(i) = COLOR_SCHEME[mod1(i, length(COLOR_SCHEME))]
 

px2(y, py, E, p)    = 2p[2] * (E - HenonHeiles.potential(0.0, y, p)) - py^2
in_section(v, E, p) = px2(v[1], v[2], E, p) > 0
pymax(y, E, p) = sqrt(max(0.0, 2 * p[2] * (E - HenonHeiles.potential(0.0, y, p))))



"Lift (y, py) to a 4D state. Offset follows sign(px) -> no phantom t=0 crossing."
function lift(v, E, p; sgn = +1)
    a = px2(v[1], v[2], E, p)
    a <= 0 && return nothing
    return [EPS_OFF, v[1], sgn * sqrt(a), v[2]]
end

"return sol, pMap, tsd"
function section_trj(v, prm)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")

    empty!(prm.yd); empty!(prm.pyd); empty!(prm.tsd)
    reinit!(prm.integ_dense, u0)
    solve!(prm.integ_dense)
    sol = prm.integ_dense.sol
    return sol, permutedims([prm.yd prm.pyd]), copy(prm.tsd)
end

"Find the smallest k with |T^k v - v| < tol.

returns named tuple (;traj, pMap, Nperiod, Tperiod)"
function minPeriodicity(v, prm; tol = 1e-8, search = 40)
    prm.nmax_dense[] = search          # dense integrator only — no interference with T

    sol, trace, ts = section_trj(v, prm)

    for i in axes(trace, 2)
        if norm(v .- trace[:, i]) < tol
            k = searchsortedfirst(sol.t, ts[i])
            return (; traj = sol.u[1:k], pMap = trace[:, 1:i],
                      Nperiod = i, Tperiod = ts[i])
        end
    end

    @warn "no closure within $search crossings" v tol
    return (; traj = sol.u, pMap = trace, Nperiod = nothing, Tperiod = nothing)
end

"""
Build the fast (residuals) and dense (trajectories) integrators.
 
Returns the nmax Refs as well -- always construct `SectionParams` from these
returned Refs, never from freshly-made ones, or writes to `prm.nmax_*` will be
invisible to the callbacks.
"""
function create_integrators(p; tmax = 20_000.0, nfast = 1, ndense = 40)
    nmax_fast, nmax_dense = Ref(nfast), Ref(ndense)
    condition(u, t, integ) = u[1]
 
    # --- fast: no trajectory saved, terminates as soon as n crossings are in ---
    yf, pyf, tsf = Float64[], Float64[], Float64[]
    affect_f!(integ) = begin
        push!(yf, integ.u[2]); push!(pyf, integ.u[4]); push!(tsf, integ.t)
        length(yf) ≥ nmax_fast[] && terminate!(integ)
    end
    integ_fast = init(ODEProblem(HenonHeiles.equations!, zeros(4), (0.0, tmax), p),
                      Vern9(); abstol = 1e-14, reltol = 1e-14,
                      save_everystep = false, save_start = false,
                      callback = ContinuousCallback(condition, affect_f!, nothing;
                                                    abstol = 1e-13))
 
    # --- dense: saves the trajectory for plotting / period detection ---
    yd, pyd, tsd = Float64[], Float64[], Float64[]
    affect_d!(integ) = begin
        push!(yd, integ.u[2]); push!(pyd, integ.u[4]); push!(tsd, integ.t)
        length(yd) ≥ nmax_dense[] && terminate!(integ)
    end
    integ_dense = init(ODEProblem(HenonHeiles.equations!, zeros(4), (0.0, tmax), p),
                       Vern9(); abstol = 1e-14, reltol = 1e-14, saveat = dt,
                       callback = ContinuousCallback(condition, affect_d!, nothing;
                                                     abstol = 1e-13))
 
    return (; integ_fast, yf, pyf, tsf, nmax_fast,
              integ_dense, yd, pyd, tsd, nmax_dense)
end


function SectionParams(E, p; tmax = 20_000.0, nfast = 1, ndense = 40)
    b = create_integrators(p; tmax, nfast, ndense)
    return SectionParams(b.integ_fast, b.yf, b.pyf, b.tsf, b.nmax_fast,
                         b.integ_dense, b.yd, b.pyd, b.tsd, b.nmax_dense,
                         E, p, tmax)
end

"n-th crossing of the section, starting from v. This is T^n(v)."
function T(v, n::Int, prm::SectionParams)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")
 
    n > prm.nmax_fast[] && (prm.nmax_fast[] = n)
 
    empty!(prm.yf); empty!(prm.pyf); empty!(prm.tsf)
    reinit!(prm.integ_fast, u0)
    solve!(prm.integ_fast)
 
    length(prm.yf) < n &&
        error("only $(length(prm.yf)) crossings in t < $(prm.tmax) (need $n)")
    return [prm.yf[n], prm.pyf[n]]
end

Fres(v, n, prm) = T(v, n, prm) - v


"Finite everywhere: TrustRegion probes outside the boundary and NaN poisons it."
function Fres_safe(v, n, prm)
    a = px2(v[1], v[2], prm.E, prm.p)
    a <= 0 && return fill(1.0 + 100 * sqrt(-a), 2)
    return Fres(v, n, prm)
end


function jacobian!(J, v, n, prm, d)
    E, p = prm.E, prm.p
    for (j, e) in enumerate(([d, 0.0], [0.0, d]))
        vp, vm = v .+ e, v .- e
        okp, okm = in_section(vp, E, p), in_section(vm, E, p)
        if okp && okm
            J[:, j] = (Fres(vp, n, prm) .- Fres(vm, n, prm)) ./ (2d)
        elseif okp
            J[:, j] = (Fres(vp, n, prm) .- Fres(v, n, prm)) ./ d       # forward
        elseif okm
            J[:, j] = (Fres(v, n, prm) .- Fres(vm, n, prm)) ./ d      # backward
        else
            error("both probes outside boundary at $v -- v is on the edge")
        end
    end
    return J
end

function jacobian(v,n,prm,d)
    J    = zeros(2, 2)
    jacobian!(J,v,n,prm,d)
end

function get_DT(v,n,prm; d=1e-7)
    I+jacobian(v,n,prm, d) 
end




function find_orbit(v0, n, prm;
                    N_max = 100, d = 1e-7, tol = 1e-11, 
                    max_backtrack = 30, dmax = 0.05, verbose = false)
 
 
    v    = collect(float.(v0))
    hist = zeros(2, N_max)
    J    = zeros(2, 2)
 
    fail(i, rn, msg) = (v = v, DT = nothing, converged = false,
                        history = hist[:, 1:max(i, 0)], resnorm = rn, comment = msg)
 
    in_section(v, prm.E, prm.p) || return fail(0, Inf, "seed outside energy boundary")
 
    r  = Fres(v, n, prm)
    rn = norm(r)
 
    for i in 1:N_max
        hist[:, i] = v
 
        if rn < tol
            jacobian!(J, v, n, prm, d)          # evaluated AT the root
            DT = J + I
            return (v = v, DT = DT, converged = true,
                    history = hist[:, 1:i], resnorm = rn,
                    comment = "|r| = $rn  det(DT) = $(det(DT))")
        end
 
        jacobian!(J, v, n, prm, d)
        κ = cond(J)
        κ > 1e12 && @warn "ill-conditioned Jacobian" iterate=i cond=κ
 
        step = J \ r
        sn   = norm(step)
        sn > dmax && (step .*= dmax / sn)     # cap: keeps Newton local
 
        λ, accepted = 1.0, false              # MUST start false
        for _ in 1:max_backtrack
            vnew = v .- λ .* step
            rnew = in_section(vnew, prm.E, prm.p) ? Fres(vnew, n, prm) : nothing
            if rnew !== nothing && norm(rnew) < rn
                v, r, rn, accepted = vnew, rnew, norm(rnew), true
                break
            end
            λ /= 2
        end
        verbose && println("  it $i  |r| = $rn  lambda = $λ  cond = $κ")
        accepted || return fail(i, rn, "line search stalled at |r| = $rn")
    end
    fail(N_max, rn, "$N_max iterations exhausted, |r| = $rn")
end


"NonlinearSolve variant. AutoFiniteDiff is mandatory: the ODE callback rejects Duals."
function solve_orbit(v0, n, prm; tol = 1e-11, maxiters = 300)
    v = collect(float.(v0))
    in_section(v, prm.E, prm.p) ||
        return (v = v, DT = nothing, converged = false, resnorm = Inf,
                history = zeros(2, 0), comment = "seed outside boundary")
 
    prob = NonlinearProblem((w, q) -> Fres_safe(w, n, q), v, prm)
    sol  = solve(prob, TrustRegion(autodiff = AutoFiniteDiff(fdtype = Val(:central)));
                 abstol = tol, maxiters)
 
    ok = SciMLBase.successful_retcode(sol) && in_section(sol.u, prm.E, prm.p)
    ok || return (v = sol.u, DT = nothing, converged = false,
                  resnorm = norm(sol.resid), history = zeros(2, 0),
                  comment = "$(sol.retcode)")
 
    DT = get_DT(sol.u, n, prm)
    return (v = sol.u, DT = DT, converged = true, resnorm = norm(sol.resid),
            history = zeros(2, 0), comment = "$(sol.retcode)  det(DT) = $(det(DT))")
end


"True if v coincides with ANY section point of an already-catalogued orbit."
function already_found(df, v, prime; tol = 1e-6)
    v === nothing && return false
    for o in eachrow(df)
        o.prime == prime || continue
        any(i -> norm([o.sec_y[i], o.sec_py[i]] .- v) < tol, eachindex(o.sec_y)) &&
            return true
    end
    return false
end


"Is the orbit invariant under S(y,py) = (y,-py)?"
function is_sym(sec; tol = 1e-6)
    Ssec = sec .* [1.0, -1.0]
    for j in axes(Ssec, 2)
        any(i -> norm(Ssec[:, j] - sec[:, i]) < tol, axes(sec, 2)) || return false
    end
    return true
end

P_py(v) = [v[1], -v[2]]


"rotate momenta and position by an angle θ"
function rotate_space(v, θ)
    c, s = cos(θ) ,sin(θ)
    return [c -s; s c]*v 
end

function rotate_state(u, θ)
    c, s = cos(θ) ,sin(θ)
    return [c*u[1]-s*u[2], s*u[1]+c*u[2], c*u[3]-s*u[4], s*u[3]+c*u[4]]
end

function is_rot_sym(sec, E, p, θ; tol = 1e-6)
    Ssec = rotate_state(lift(sec, E, p),θ)
    for j in axes(Ssec, 2)
        any(i -> norm(Ssec[:, j] - sec[:, i]) < tol, axes(sec, 2)) || return false
    end
    return true
end


function rotated_root(v, prm, θ; search = 8, margin = 1e-8)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && return nothing

    old = prm.nmax_dense[]
    prm.nmax_dense[] = search          # a few crossings, so we can be choosy
    try
        empty!(prm.yd); empty!(prm.pyd); empty!(prm.tsd)
        reinit!(prm.integ_dense, rotate_state(u0, θ))
        solve!(prm.integ_dense)
        for i in eachindex(prm.yd)
            w = [prm.yd[i], prm.pyd[i]]
            px2(w[1], w[2], prm.E, prm.p) > margin && return w
        end
        return nothing                 # all crossings degenerate
    finally
        prm.nmax_dense[] = old
    end
end

"""
Extend the catalogue using the symmetry group -- no Newton, one integration each.
 
S costs nothing at all; each C_3 rotation costs a single flow to the section.
Both are orders of magnitude cheaper than a fresh solve.
"""
function symmetrise_E!(df, prm; tol = 1e-6, verify = true)
    base = collect(eachrow(df))
    nid  = isempty(df) ? 0 : maximum(df.id)
 
    for o in base
        # --- time reversal: free ---
        if !is_sym(sec_matrix(o); tol)
            v1 = P_py([o.y, o.py])
            if !already_found(df, v1, o.prime; tol)
                verify && @printf("  S      id=%d  |F(v)| = %.2e\n",
                                  o.id, norm(Fres(v1, o.prime, prm)))
                new = analyse_seed(v1, o.prime, prm; id = nid + 1, origin = "S(id$(o.id))")
                if new !== nothing
                    nid += 1
                    push!(df, merge(new, (; id = nid)))
                end
            end
        end
 
        # --- C_3 rotations: one integration each ---
        for k in 1:2
            v1 = try
                rotated_root([o.y, o.py], prm, k * 2π / 3)
            catch e
                @warn "rotation failed" id = o.id k exception = e
                continue
            end
            v1 === nothing && continue          # <-- must precede already_found
            already_found(df, v1, o.prime; tol) && continue
            verify && @printf("  C3^%d   id=%d  |F(v)| = %.2e\n",
                              k, o.id, residual_norm(v1, o.prime, prm))
            new = analyse_seed(v1, o.prime, prm;
                               id = nid + 1, origin = "C3^$k(id$(o.id))")
            if new !== nothing
                nid += 1
                push!(df, merge(new, (; id = nid)))
            end
        end
    end
    return df
end


function symmetrise(df, p; tmax = 10_000.0, ndense = 40, tol = 1e-6)
    out = orbit_table()
    for g in groupby(df, :E)
        sub = DataFrame(g)
        prm_E = SectionParams(first(sub.E), p; tmax, ndense = ndense)
        symmetrise_E!(sub, prm_E; verify = false)
        append!(out, sub)
    end
    return out
end

"cleans the data to contain no repetitions"
function dedup(df; tol = 1e-6)
    keep = trues(nrow(df))
    for i in 2:nrow(df)
        for j in 1:i-1
            keep[j] || continue
            df.prime[i] == df.prime[j] || continue
            same = any(eachindex(df.sec_y[i])) do a
                any(eachindex(df.sec_y[j])) do b
                    hypot(df.sec_y[i][a] - df.sec_y[j][b],
                          df.sec_py[i][a] - df.sec_py[j][b]) < tol
                end
            end
            if same
                keep[i] = false
                break
            end
        end
    end
    return df[keep, :]
end

"""
    dedup_all(df; tol = 1e-6)

Collapse repeated detections of the same orbit. Orbits are compared only within
the same energy — section points at different E can coincide numerically without
being the same orbit.

Within each energy group, rows are sorted by `prime` so the shortest orbit is
kept as the representative, and its repeats at higher `n` are dropped.
`nhits` records how many raw rows collapsed into each entry (a proxy for basin size).
"""
function dedup_all(df; tol = 1e-6)
    isempty(df) && return copy(df)
    groups = DataFrame[]
    for g in groupby(df, :E)
        push!(groups, dedup(sort(DataFrame(g), [:prime, :resnorm]); tol = tol))
    end
    sort!(vcat(groups...), [:E, :prime, :T])
end

"creating the initial points from where the search starts from"
function section_grid(E, p; ny = 10, npy = 10, margin = 0.03)
    roots = real.(filter(r -> abs(imag(r)) < 1e-10, HenonHeiles.limit_of_initial_y0(E, p)))
    sort!(roots)
    ymin, ymax = roots[1], roots[2]      # the two turning points bounding the well
    dy = margin * (ymax - ymin)
    seeds = Vector{Float64}[]
    for y in range(ymin + dy, ymax - dy, length = ny)
        pm = pymax(y, E, p)
        pm <= 0 && continue
        for s in range(0 + margin, 1 - margin, length = npy)
            push!(seeds, [y, s * pm])
        end
    end
    seeds
end


"returns a row to the orbit"
function analyse_seed(v0, n, prm; id = 0, origin = "seed",
                      tol = 1e-11, search = 40)
    res = solve_orbit(v0, n, prm; tol)
    res.converged || return nothing
 
    orb = minPeriodicity(res.v, prm; search)
    orb.Nperiod === nothing && return nothing
 
    DT    = get_DT(res.v, orb.Nperiod, prm)
    τ     = tr(DT)
    index = kind_index(τ)
 
    return (; E = prm.E, n, seed_y = v0[1], seed_py = v0[2],
              y = res.v[1], py = res.v[2],
              prime = orb.Nperiod, T = orb.Tperiod,
              trace = τ, detDT = det(DT),
              resnorm = res.resnorm, iters = size(res.history, 2),
              sec_y = orb.pMap[1, :], sec_py = orb.pMap[2, :],
              history_y = res.history[1, :], history_py = res.history[2, :],
              traj_x = [u[1] for u in orb.traj], traj_y = [u[2] for u in orb.traj],
              class = KIND_LABEL[index], index, id, origin)
end


"""
    Takes 
"""
function sweep!(df::AbstractDataFrame, seeds, n, prm; tol = 1e-11, verbose = false)
    nid = isempty(df) ? 0 : maximum(df.id)
    @showprogress dt = 1 desc = "E=$(round(prm.E; digits = 4)) n=$n" for v0 in seeds
        orb = try
            analyse_seed(v0, n, prm; id = nid + 1, tol)
        catch e
            verbose && @warn "seed failed" seed = v0 exception = e
            continue
        end
        orb === nothing              && continue
        orb.prime == n               || continue
        nid += 1
        push!(df, merge(orb, (; id = nid)))
    end
    return df
end


"""
Sweep all periods 1:nmax at one energy. Serial by construction: the n loop
shares a df so that `already_found` can reject rediscoveries, and low-n orbits
must be catalogued before high-n so the divisor filter has something to match.

returns (;df, timing)
"""
function sweep_energy(E, ns, p; tmax = 5000.0, ny = 5, npy = 5,
                      ndense = 40, seeds = nothing, symmetrise = false)
    prm = SectionParams(E, p; tmax, nfast = 1, ndense)
    seeds === nothing && (seeds = section_grid(E, p; ny, npy))
    df  = orbit_table()
    timing = DataFrame(E = Float64[], n = Int[], nseeds = Int[],
                       found = Int[], secs = Float64[])

    for n in ns
        before = nrow(df)
        secs   = @elapsed sweep!(df, seeds, n, prm)
        push!(timing, (E, n, length(seeds), nrow(df) - before, secs))
    end

    symmetrise && symmetrise_E!(df, prm; verify = false)
    return (; df, timing)
end

"""
Run `sweep_energy` across energies on all available threads.

One SectionParams per task, created inside the task -- never shared, because
the integrators and their callback buffers are mutable state.
"""
function run_sweep_mt(Es, ns, p; tmax = 5000.0, ny = 5, npy = 5,
                      ndense = 40, symmetrise = true)
    println("running sweep, parallel computation\nSymmitrise=$symmetrise\n")

    tasks = map(Es) do E
        Threads.@spawn sweep_energy(E, ns, p; tmax, ny, npy, ndense,symmetrise=symmetrise)
    end
    results = fetch.(tasks)


    df     = reduce(vcat, r.df     for r in results)
    timing = reduce(vcat, r.timing for r in results)


    sort!(df, [:E, :prime, :T])
    return (; df, timing)
end


"""
Trace periodic orbtis for various Energies Es.
This function starts at a given max(Es). Finds many periodic orbits.
Then uses the found periodic orbits for the next energy to find it again.

returns (;df, timing)
"""
function follow_orbits(Es, ns, p; tmax = 5000.0, ny_init = 15, npy_init = 10,
                      ndense = 40, symmetrise = true)
    sort!(Es; rev=true)
    println("Following periodic orbits from E=$(Es[1]) to $(Es[end])\n")

    initial_sweep   =       sweep_energy(Es[1], ns, p; tmax, ny=ny_init, npy=npy_init, ndense, symmetrise)
    red_init        =       dedup_all(initial_sweep.df; tol = 1e-6)
    seeds = [[red_init.y[i], red_init.py[i]] for i in 1:nrow(red_init)]    

    all_df          =       [initial_sweep.df]
    all_timing      =       [initial_sweep.timing]
    for e in Es[2:end]
        r = sweep_energy(e, ns, p; tmax, ndense, seeds, symmetrise)
        push!(all_df, r.df)
        push!(all_timing, r.timing)
        reduced = dedup_all(r.df)
        seeds = [[reduced.y[i], reduced.py[i]] for i in 1:nrow(reduced)]    
    end
    df = reduce(vcat, all_df)
    timing =reduce(vcat,all_timing)


    sort!(df, [:E, :prime, :T])
    sort!(timing, [:E])
    return (;df, timing)
end




function report(df)
    println("="^78)
    @printf("%3s %5s %6s %9s %9s %11s %11s  %-20s %s\n",
            "id", "n", "prime", "y", "py", "T", "tr(DT)", "class", "origin")
    println("-"^78)
    for o in eachrow(df)
        @printf("%3d %5d %6d %9.5f %9.5f %11.4f %11.4f  %-20s %s\n",
                o.id, o.n, o.prime, o.y, o.py, o.T, o.trace, o.class, o.origin)
    end
    println("="^78)
    bad = filter(o -> abs(o.detDT - 1) > 1e-4, eachrow(df))
    isempty(bad) || @warn "rows with det(DT) far from 1 -- check the FD step" ids = [o.id for o in bad]
end
 


function few_orbits()
    symmetrise   = false
    p            = (1.0, 1.0, 1.0)
    E            = 0.1144
    n            = 6              # crossings the dense integrator may take
    ndense       = 10
    tmax         = 10_000.0
    
    prm = SectionParams(E, p; tmax, nfast = n, ndense)
    seeds = section_grid(E, p; ny=6, npy=3)
    df  = orbit_table()
    timing = DataFrame(E = Float64[], n = Int[], nseeds = Int[],
                       found = Int[], secs = Float64[])

    before = nrow(df)
    secs   = @elapsed sweep!(df, seeds, n, prm)
    push!(timing, (E, n, length(seeds), nrow(df) - before, secs))

    symmetrise && symmetrise_E!(df, prm; verify = false)
    return (; df, timing)
end

orbit_table()

function main()
    res1  = few_orbits()
    orbti, _ =res.df, res.timing

    for o in eachrow(orbits)
        println()
    
    end

end

res = main()