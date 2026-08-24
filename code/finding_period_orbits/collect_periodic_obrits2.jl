import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2,
      NonlinearSolve, ADTypes, DataFrames, Dates, ProgressMeter, Printf, ColorSchemes
include("../styles/makie_theme.jl")
set_theme!(QUARTO_THEME)
include("../models/henon_heiles.jl")
using .HenonHeiles

const CONFIG_DIR = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
const DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
const FIG_DIR    = joinpath(@__DIR__, "../../figures/henon-heiles/periodic-orbits/")
data_file        = joinpath(DATA_DIR, "E0.1127-T10000.0-py0.0-n256.jld2")

configurations = JSON3.read(read(CONFIG_DIR, String))
cfg            = configurations.explore
param          = [Float64(cfg.a.value), Float64(cfg.m.value), Float64(cfg.w.value)]
dt             = Float64(cfg.dt.value)
x0             = Float64(cfg.x0.value)


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
 

# --- background section ---
data = load(data_file, "results")
y_all, py_all = Float32[], Float32[]
for d in data
    append!(y_all,  d.sec_y)
    append!(py_all, d.sec_py)
end

const EPS_OFF = 1e-15

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
P_pypx(u) = [u[1], u[2], -u[3], -u[4]]
P_x(u) = [-u[1], u[2], -u[3], u[4]]


"rotate momenta and position by an angle θ"
function rotate_space(v, θ)
    c, s = cos(θ) ,sin(θ)
    return [c -s; s c]*v 
end

function rotate_state(u, θ)
    c, s = cos(θ) ,sin(θ)
    return [c*u[1]-s*u[2], s*u[1]+c*u[2], c*u[3]-s*u[4], s*u[3]+c*u[4]]
end

function is_rot_sym(sec, E, p; tol = 1e-6)
    Ssec = rotate_state(lift(sec, E, p))
    for j in axes(Ssec, 2)
        any(i -> norm(Ssec[:, j] - sec[:, i]) < tol, axes(sec, 2)) || return false
    end
    return true
end

function rotated_root(v, prm, θ)
    u_rot = rotate_state(lift(v, prm.E, prm.p), θ)
    empty!(prm.yd); empty!(prm.pyd); empty!(prm.tsd)
    reinit!(prm.integ_dense, u_rot)
    solve!(prm.integ_dense)
    isempty(prm.yd) && error("rotated state never crossed the section")
    return [prm.yd[1], prm.pyd[1]]
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


" 
E::Float64; n::Int
seed_y::Float64; seed_py::Float64
y::Float64; py::Float64
prime::Int; T::Float64
trace::Float64; detDT::Float64
resnorm::Float64; iters::Int
sec_y::Vector{Float64}; sec_py::Vector{Float64}
history_y::Vector{Float64}; history_py::Vector{Float64}
traj_x::Vector{Float64}; traj_y::Vector{Float64}
"
orbit_table() = DataFrame(Row[])

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



function sweep!(df, seeds, n, prm; tol = 1e-11, verbose = false)
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


function run_sweep(Es, ns, p; tmax = 5000.0, ny = 5, npy = 5)
    df = orbit_table()
    timing = DataFrame(E = Float64[], n = Int[], nseeds = Int[],
                       found = Int[], secs = Float64[])
    for E in Es, n in ns
        bff = create_integrators(p; tmax, nmax = Ref(n))
        prm = (; integ = bff.integ, y = bff.y, py = bff.py, ts = bff.ts,
                 E, p, n, tmax)
        seeds = section_grid(E, p; ny, npy)
        before = nrow(df)  
        secs = @elapsed sweep!(df, seeds, prm)
        push!(timing, (E, n, length(seeds), nrow(df) - before, secs))
    end
    (; df, timing)
end


"cleans the data to contain no repetitions"
function dedup(df; tol = 1e-6)
    keep = trues(nrow(df))
    sec = permutedims([df.sec_y df.sec_py])
    for i in 2:nrow(df)
        for j in 1:i-1
            keep[j] || continue
            
            if any(norm(a .- b) < tol for a in sec[i], b in sec[j])
                keep[i] = false
                break
            end
        end
    end
    df[keep, :]
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
 
# =====================================================================
#  Plotting
# =====================================================================
 
function plot_config(df, prm; title = "configuration space")
    r      = range(-1.0, 1.0, length = 200)
    epot   = [HenonHeiles.potential(x, y, prm.p) for x in r, y in r]
    levels = logrange(5.0 * 0.009, 6.9 * 0.089, 7)
 
    fig = Figure(size = (1400, 900))
    ax  = Axis(fig[1, 1], xlabel = "x", ylabel = "y",
               title = title, aspect = DataAspect())
    contour!(ax, r, r, epot; levels, colormap = :hsv, labels = true)
    contour!(ax, r, r, epot; levels = [prm.E], color = C_CREAM, linewidth = 2)
 
    for o in eachrow(df)
        lines!(ax, o.traj_x, o.traj_y; color = pick_color(o.id),
               linestyle = KIND_LS[o.index],
               label = "id$(o.id) n=$(o.prime) T=$(round(o.T; digits = 1))")
    end
    kinds = sort(unique(df.index))
    elems = [LineElement(linestyle = KIND_LS[k], linewidth = 2) for k in kinds]
    Legend(fig[2, 1], elems, KIND_LABEL[kinds], "stability";
           orientation = :horizontal, framevisible = false)    
    return fig, ax
end
 
function plot_section(df, prm; seeds = nothing, title = "surface of section")
    fig = Figure(size = (1400, 900))
    ax  = Axis(fig[1, 1], xlabel = "y", ylabel = "p_y", title = title)
 
    scatter!(ax, y_all, py_all; color = (:grey, 0.5), markersize = 1.5)
    y_max, py_max = HenonHeiles.section_boundary_ranges(prm.E, prm.p, 120)
    scatter!(ax, HenonHeiles.section_boundary(y_max, py_max);
             color = C_CREAM, markersize = 4)
 
    seeds === nothing || scatter!(ax, first.(seeds), last.(seeds);
                                  color = (:white, 0.35), markersize = 5)
 
    for o in eachrow(df)
        scatter!(ax, o.sec_y, o.sec_py; color = pick_color(o.id), markersize = 12,
                 marker = KIND_MS[o.index],
                 label = "id$(o.id) n=$(o.prime) $(o.class)")
    end
    kinds = sort(unique(df.index))
    elems = [MarkerElement(marker = KIND_MS[k], markersize = 12, color = :white) for k in kinds]
    Legend(fig[2, 1], elems, KIND_LABEL[kinds], "stability";
           orientation = :horizontal, framevisible = false)
    return fig, ax
end

function plot_orbits(df, prm; seeds = nothing)
    fig = Figure(size = (1800, 950))

    ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = "y",
               title = "configuration space", aspect = DataAspect())
    ax2 = Axis(fig[1, 2], xlabel = "y", ylabel = "p_y",
               title = "surface of section")

    # --- backdrop ---
    r    = range(-1.0, 1.0, length = 200)
    epot = [HenonHeiles.potential(x, y, prm.p) for x in r, y in r]
    contour!(ax1, r, r, epot; levels = [prm.E], color = C_CREAM, linewidth = 2)
    scatter!(ax2, y_all, py_all; color = (:grey, 0.4), markersize = 1.5)
    y_max, py_max = HenonHeiles.section_boundary_ranges(prm.E, prm.p, 120)
    scatter!(ax2, HenonHeiles.section_boundary(y_max, py_max);
             color = C_CREAM, markersize = 3)
    seeds === nothing || scatter!(ax2, first.(seeds), last.(seeds);
                                  color = (:white, 0.25), markersize = 4)

    # --- shared colour scale on T ---
    Ts   = df.T
    crange = (minimum(Ts), maximum(Ts))
    cmap   = :viridis

    for o in eachrow(df)
        c = get(colorschemes[cmap], (o.T - crange[1]) / (crange[2] - crange[1]))
        lines!(ax1, o.traj_x, o.traj_y; color = c, linestyle = KIND_LS[o.index],
               linewidth = 2)
        scatter!(ax2, o.sec_y, o.sec_py; color = c, marker = KIND_MS[o.index],
                 markersize = 14, strokewidth = 0.5, strokecolor = :black)
    end

    Colorbar(fig[1, 3]; limits = crange, colormap = cmap, label = "period T")

    # --- legend for the categorical channel only ---
    kinds = sort(unique(df.index))
    elems = [[LineElement(linestyle = KIND_LS[k], linewidth = 2),
              MarkerElement(marker = KIND_MS[k], markersize = 12, color = :white)]
             for k in kinds]
    Legend(fig[2, 1:3], elems, KIND_LABEL[kinds], "stability";
           orientation = :horizontal, framevisible = false)

    return fig, ax1, ax2
end


p            = (1.0, 1.0, 1.0)
E0           = 0.1127
nmax         = 6              # highest period to search (2D Newton degrades past ~6)
nmax_search  = 8              # crossings the dense integrator may take
tmax         = 20_000.0
ny, npy      = 10, 5
 
prm    = SectionParams(E0, p; tmax, nfast = 1, ndense = nmax_search)
seeds  = section_grid(E0, p; ny, npy)
orbits = orbit_table()
 
println("seeds: $(length(seeds))   periods: 1:$nmax")
for n in 1:nmax
    sweep!(orbits, seeds, n, prm)
end


# f = joinpath(SAVE_DATA_DIR, "orbits_E0.01-0.1127_20260816-0742.jld2")   

res = dedup_all(orbits)

base = collect(eachrow(res))
# for i in 1:2, orb in base
#     v1 = try
#         rotated_root([orb.y, orb.py], prm, i*2π/3)
#     catch e
#         @warn "rotation failed" seed=(orb.y, orb.py) exception=e
#         continue
#     end
#     @show i*2/3, norm(Fres(v1, orb.prime, prm)), orb.id
#     new = analyse_seed(v1, 1, prm; id= orb.id)
#     # new.id=orb.id
#     new === nothing || push!(orbits, new)
# end

@showprogress dt = 1 desc = "symmertrizing" for orb in base
    sec = permutedims([orb.sec_y orb.sec_py])     # rebuild the 2 x k matrix
    is_sym(sec; tol = 1e-6) && continue            # already closed under S — skip
    v1  = [orb.y, -orb.py]
    new = analyse_seed(v1, orb.prime, prm; id = orb.id)
    new === nothing || push!(orbits, new)
end
res = dedup_all(orbits)


pltorb, axorb = plot_orbits(orbits,prm)
display(pltorb)
# cs   = [C_RED, C_PURPLE, C_CREAM, C_TEAL, C_ORANGE, C_GOLD,C_GREEN, C_BLUE,C_PINK,C_GREY]

# r      = range(-1.0, 1.0, length=120)
# levels = logrange(5.0*0.009, 6.9*0.089, 7)
# epot   = [HenonHeiles.potential(x,y, param) for x in r, y in r]

# fig_conf = Figure(size = (1400, 900))
# ax_conf  = Axis(fig_conf[1, 1], xlabel = "x", ylabel = "y",
#              title = "config space, n = 1", aspect = DataAspect())
# contour!(ax_conf, r, r, epot, labels=true, levels=levels, colormap=:hsv, colorscale=identity)
# for (j,o) in enumerate(eachrow(orbits))
#     # j >4 && lines!(ax_conf, o.traj_x, o.traj_y, color = color_scheme[o.id],
#     #        label = "T=$(round(o.T; digits=1)), $(o.class)", linestyle=KIND_LS[o.index])
#     # j==4 && lines!(ax_conf, o.traj_x, o.traj_y, color = color_scheme[o.id],
#     #        label = "T=$(round(o.T; digits=1)), $(o.class)", linestyle=(:dash))
#     o.primes==1 && lines!(ax_conf, o.traj_x, o.traj_y, color = cs[o.prime],
#             label = "T=$(round(o.T; digits=1)), $(o.class)", linestyle=KIND_LS[o.index])
# end #round(x; digits = 3)
# Legend(fig_conf[2, 1], ax_conf; orientation = :horizontal, framevisible = false)
# display(GLMakie.Screen(), fig_conf)

# fig_p = Figure(size=(1400,900))
# ax_p = Axis(fig_p[1, 1], xlabel = "y", ylabel = "p_y",
#                title = "section, n = 1, prime = 1")
    

# y_max, py_max = HenonHeiles.section_boundary_ranges(E0, param, 120)
# boundary      = HenonHeiles.section_boundary(y_max, py_max) 
# scatter!(ax_p, y_all, py_all, color = (:grey, 0.5), markersize = 1.5)
# scatter!(ax_p, boundary, color = C_CREAM, markersize = 4)
# # scatterlines!(ax_p, res.history[1, :], res.history[2, :],
#                     # color = C_GOLD, markersize = 8, label = "Newton path")
# #scatter!(ax_p, first.(initial_search), last.(initial_search), color=C_TEAL, markersize=8)
# for (j, o) in enumerate(eachrow(orbits))
#     # j>4 && scatter!(ax_p, o.sec_y, o.sec_py,
#     #          color = color_scheme[o.id], markersize = 12,
#     #          marker = KIND_MS[o.index],
#     #          label = KIND_LABEL[o.index])
#     scatter!(ax_p, o.sec_y, o.sec_py,
#              color = cs[o.prime], markersize = 12,
#              marker = KIND_MS[o.index],
#              label = KIND_LABEL[o.index])
# end
# Legend(fig_p[2, 1],
#        [MarkerElement(marker = m, color =C_CREAM) for m in KIND_MS],
#        KIND_LABEL; orientation = :horizontal, framevisible = false)
# Legend(fig_p[1, 2],
#        [MarkerElement(marker = m, color =cs) for m in KIND_MS],
#        KIND_LABEL; orientation = :horizontal, framevisible = false)
# display(GLMakie.Screen(), fig_p)




# sol=Orbit_finder(v0, n, E0; theta=2/3)

# save(joinpath(FIG_DIR, "phsp-two-orbits-by-Sreflection-AndBase-$(E0)-n$nmax.png"),fig_p; px_per_unit = 2)
# save(joinpath(FIG_DIR, "conf-two-orbits-by-Sreflection-AndBase-$(E0)-n$nmax.png"), fig_conf; px_per_unit = 2)

# println("\n")
# for (j,o) in enumerate(eachrow(orbits))
#     j <= 4 && println("$j) Base orbit: $(o.id)\n y=$(o.y), py=$(o.py)\n|r|=$(o.resnorm)")
# end
