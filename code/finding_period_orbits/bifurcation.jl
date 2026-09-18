import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2,
      NonlinearSolve, ADTypes, DataFrames, Dates, ProgressMeter, Printf, ColorSchemes
include("../styles/makie_theme.jl")
include("../models/henon_heiles.jl")
using .HenonHeiles

BLAS.set_num_threads(1)          # your linear algebra is 2x2; BLAS threads only compete

const CONFIG_DIR = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
const DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
const FIG_DIR    = joinpath(@__DIR__, "../../figures/henon-heiles/sep16-latex/")
const SAVE_DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/bifurcation/")



const RGRID  = range(-2.0, 2.0, length = 840)
const EPOT   = [HenonHeiles.potential(x, y, (1.0,1.0,1.0)) for x in RGRID, y in RGRID]
const LEVELS = [0.001, 0.01, 0.05,  0.1, 0.14,  0.166666,0.24,  0.6, 1.0, 1.5,2.3]# collect(logrange(5.0*0.009, 6.9*0.089, 10))


if false
    save_name= "pot_contour.png"
    save_as=joinpath(FIG_DIR, save_name)
    set_style!(:print)
    Fpot = Figure(size=(1300, 900))
    axpot = Axis(Fpot[1,1], ylabel="y", xlabel="x")
    contour!(axpot, RGRID, RGRID, EPOT, colormap=:viridis, labels=true,levels=LEVELS, linewidth=2.5)
    xlims!(-1.25,1.25)
    ylims!(-1.25,1.25)
    save(save_as,Fpot, px_per_unit = 2)
    display(Fpot)
end


const CC_TOL  = 1e-13
const INT_TOL = 1e-14
const PMAP_ROOT_TOL = 1e-11
const PMAP_PRIME_TOL = 1e-10
const OFFSET_X_SURFACESEC = 1e-13

mutable struct SectionParams1{IF, P}
    integ_fast  :: IF
    yf          :: Vector{Float64}
    pyf         :: Vector{Float64}
    tsf         :: Vector{Float64}
    nmax_fast   :: Ref{Int}
 
    E           :: Float64
    p           :: P
    tmax        :: Float64
end

mutable struct SectionParams{IF, ID, P}
    integ_fast  :: IF
    yf          :: Vector{Float64}
    pyf         :: Vector{Float64}
    tsf         :: Vector{Float64}
    nmax_fast   :: Ref{Int}
 
    integ_dense :: ID
    yd          :: Vector{Float64}
    pyd         :: Vector{Float64}
    tsd         :: Vector{Float64}
    nmax_dense  :: Ref{Int}

    E           :: Float64
    p           :: P
    tmax        :: Float64
end

set_energy!(prm::SectionParams, E) =(prm.E = E) 



const Row = @NamedTuple begin
    E::Float64 
    v::Vector{Float64}
    T::Float64
    str:: String
end
 
""" 
    E::Float64
    v::Vector{Float64}
    T::Float64
    str:: String
"""
orbit_table() = DataFrame(Row[])




"returns -> 1 := elliptic, 2 := hyperbollic (stable), 3 := hyperbollic (unstable), 4 := parabolic,"
function kind_index(τ; ε = 1e-6)
    abs(abs(τ) - 2) < ε && return 4   # parabolic
    abs(τ) < 2 && return 1            # elliptic
    τ > 0 ? 2 : 3                     # hyperbollic (unstable, stable)
end

const KIND_LS    = [:solid, :dash, :dashdot, :dot]
const KIND_MS    = [:circle, :xcross, :diamond, :utriangle]
const KIND_LABEL = ["elliptic", "hyperbolic", "inverse hyperbolic", "parabolic"]




px2(y, py, E, p)    = 2p[2] * (E - HenonHeiles.potential(0.0, y, p)) - py^2
in_section(v, E, p) = px2(v[1], v[2], E, p) > 0
pymax(y, E, p) = sqrt(max(0.0, 2 * p[2] * (E - HenonHeiles.potential(0.0, y, p))))



boundary(E,p) = HenonHeiles.section_boundary(HenonHeiles.section_boundary_ranges(E, p, 100)...)

"Lift (y, py) to a 4D state. Offset follows sign(px) -> no phantom t=0 crossing."
function lift(v, E, p; sgn = +1)
    a = px2(v[1], v[2], E, p)
    a <= 0 && return nothing
    return [OFFSET_X_SURFACESEC, v[1], sgn * sqrt(a), v[2]]
end





function get_traj(u0, t;p=(1.0,1.0,1.0), abstol=INT_TOL, reltol=INT_TOL)
    prob = ODEProblem(HenonHeiles.equations!, u0, (0.0, t), p)
    sol  = solve(prob, Vern9(); abstol=abstol, reltol=reltol)
    return sol.u
end

"""
    flow ϕₜ takes u(0) to u(t)
The equations are defined in HenonHeiles.equations!
p = (1,1,1) i.e. 
    α = 1
    m = 1
    w = 1
"""
function flow(u0, t; p=(1.0,1.0,1.0), abstol = INT_TOL, reltol = INT_TOL)
    return get_traj(u0, t; p=p, abstol = abstol, reltol = reltol)[end]
end

function monodromy(u0, t; p=(1.0,1.0,1.0), d=1e-7, abstol = INT_TOL, reltol = INT_TOL)

    M   = zeros(4, 4)
    for j in 1:4
        e = zeros(4); e[j] = d
        M[:, j] = (flow(u0 .+ e, t; p=p, abstol=abstol, reltol=reltol) .-
                   flow(u0 .- e, t; p=p, abstol=abstol, reltol=reltol)) ./ (2d)
    end
    return M
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
function minPeriodicity(v, prm; pmap_prime_tol = PMAP_PRIME_TOL, search = 40)
    prm.nmax_dense[] = search          # dense integrator only — no interference with T

    sol, trace, ts = section_trj(v, prm)

    for i in axes(trace, 2)
        if norm(v .- trace[:, i]) < pmap_prime_tol
            k = searchsortedfirst(sol.t, ts[i])
            return (; traj = sol.u[1:k], pMap = trace[:, 1:i],
                      Nperiod = i, Tperiod = ts[i])
        end
    end

    @warn "no closure within $search crossings" v pmap_prime_tol
    return (; traj = sol.u, pMap = trace, Nperiod = nothing, Tperiod = nothing)
end

"""
Build the fast (residuals) and dense (trajectories) integrators.
 
Returns the nmax Refs as well -- always construct `SectionParams` from these
returned Refs, never from freshly-made ones, or writes to `prm.nmax_*` will be
invisible to the callbacks.
"""
function integrator(p; 
    tmax = 20_000.0, kind=:fast, n = 1, cc_tol = CC_TOL, int_tol=INT_TOL, save_everystep = false, save_start = false)
    nmax = Ref(n)
    condition(u, t, integ) = u[1]
 
    # --- fast: no trajectory saved, terminates as soon as n crossings are in ---
    y, py, ts = Float64[], Float64[], Float64[]
    affect_f!(integ) = begin
        push!(y, integ.u[2]); push!(py, integ.u[4]); push!(ts, integ.t)
        length(y) ≥ n[] && terminate!(integ)
    end
    cont_callback   = ContinuousCallback(condition, affect_f!, nothing;
                                                    abstol = cc_tol)
    prob            = ODEProblem(HenonHeiles.equations!, zeros(4), (0.0, tmax), p)
    if kind == :dense
        integ      = init(prob, Vern9(); 
                            abstol = 1e-14, reltol = int_tol,
                            save_everystep = save_everystep, save_start = save_start,
                            callback = cont_callback)
    elseif kind == :fast
        integ      = init(prob, Vern9();
                            abstol = 1e-14, reltol = int_tol,
                            save_everystep = false, save_start = false,
                            callback = cont_callback)
    else
        error("kind must be :fast or :dense")
    end
    return (; integ, y, py, ts, nmax)
end

"""
Build the fast (residuals) and dense (trajectories) integrators.
 
If you want to look at the orbits in configuration space, set save_everystep = true in the dense integrator, and then use `section_trj` to get the trajectory.
"""
function create_integrators(p; 
             tmax = 20_000.0, nfast = 1, ndense = 40, 
             cc_tol = CC_TOL, int_tol = INT_TOL, 
             save_everystep = save_everystep, save_start = save_start)

    condition(u, t, integ) = u[1]
 
    # --- fast: no trajectory saved, terminates as soon as n crossings are in ---
    i_fast = integrator(p; tmax, kind=:fast, n=nfast, cc_tol, int_tol, save_everystep = false, save_start = false)
 
    # --- dense: saves the trajectory for plotting / period detection ---
    i_dense = integrator(p; tmax, kind=:dense, n=ndense, cc_tol, int_tol, save_everystep = save_everystep, save_start = save_start)
 
    return (; integ_fast = i_fast.integ, yf = i_fast.y, pyf = i_fast.py, tsf = i_fast.ts, nmax_fast = i_fast.nmax, integ_dense = i_dense.integ, yd = i_dense.y, pyd = i_dense.py, tsd = i_dense.ts, nmax_dense = i_dense.nmax)
end


function SectionParams(E, p, num_int; 
          tmax = 20_000.0, nfast = 1, ndense=40, 
          cc_tol = CC_TOL, int_tol = INT_TOL, save_everystep = true, save_start = true)
    # if num_int == :fast
    #     b = root_integrator(p; tmax, nfast, cc_tol = cc_tol, int_tol = int_tol, save_everystep = false, save_start = false)
    #     return SectionParams1(b.integ_fast, b.yf, b.pyf, b.tsf, b.nmax_fast, E, p, tmax)
    # elseif num_int == :dense
    #     b = root_integrator(p; tmax, nfast, cc_tol = cc_tol, int_tol = int_tol, save_everystep = seve_everystep, save_start = save_start)
    #     return SectionParams1(b.integ_fast, b.yf, b.pyf, b.tsf, b.nmax_fast, E, p, tmax)
    if num_int == :both
        b = create_integrators(p; tmax, nfast=nfast, ndense=ndense, cc_tol, int_tol, save_everystep = save_everystep, save_start = save_start)
        return SectionParams(b.integ_fast, b.yf, b.pyf, b.tsf, b.nmax_fast, b.integ_dense, b.yd, b.pyd, b.tsd, b.nmax_dense, E, p, tmax)
    else
        error("num_int must be :both")
    end
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
function Fres_safe(v, n, prm::SectionParams)
    a = px2(v[1], v[2], prm.E, prm.p)
    a <= 0 && return fill(1.0 + 100 * sqrt(-a), 2)
    return Fres(v, n, prm)
end


function jacobian!(J, v, n, prm::SectionParams, d)
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
                    N_max = 100, d = 1e-7, tol = PMAP_ROOT_TOL, 
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
            return (v = v, converged = true,
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


"""
    return (v = sol.u, DT = nothing, converged = true, resnorm = norm(sol.resid),
            history = zeros(2, 0), comment = (sol.retcode)  det(DT) = (det(DT)))

NonlinearSolve variant. AutoFiniteDiff is mandatory: the ODE callback rejects Duals.
"""
function solve_orbit(v0, n, prm; tol = PMAP_ROOT_TOL, maxiters = 300)
    v = collect(float.(v0))
    in_section(v, prm.E, prm.p) ||
        return (v = v, converged = false, resnorm = Inf,
                history = zeros(2, 0), comment = "seed outside boundary")
 
    prob = NonlinearProblem((w, q) -> Fres_safe(w, n, q), v, prm)
    sol  = solve(prob, TrustRegion(autodiff = AutoFiniteDiff(fdtype = Val(:central)));
                 abstol = tol, maxiters)
 
    ok = SciMLBase.successful_retcode(sol) && in_section(sol.u, prm.E, prm.p)
    ok || return (v = sol.u, converged = false,
                  resnorm = norm(sol.resid), history = zeros(2, 0),
                  comment = "$(sol.retcode)")
 
    # DT = get_DT(sol.u, n, prm)
    return (v = sol.u, converged = true, resnorm = norm(sol.resid),
            history = zeros(2, 0), comment = "$(sol.retcode)")
end

"""
    return (; E = prm.E, v = res.v, str, T)

returns a row to the orbit
"""
function analyse_seed(v0, n, prm; str = "seed",
                      pmap_root_tol = PMAP_ROOT_TOL, pmap_prime_tol= PMAP_PRIME_TOL, maxiters=300)
    res = solve_orbit(v0, n, prm, tol=pmap_root_tol, maxiters=maxiters)
    # println(res.resnorm, res.comment)
    res.converged || return nothing
    
    min_period = minPeriodicity(res.v, prm; pmap_prime_tol = pmap_prime_tol, search = 40)
    min_period.Nperiod === nothing && return nothing
    T = min_period.Tperiod
    return (; E = prm.E, v = res.v, str, T)
end


function get_eigenvals(M)
    c = tr(M)/2 - 1
    δ = sqrt(Complex(c^2 - 1))
    (c + δ, c - δ)
end


"True if v coincides with ANY section point of an already-catalogued orbit."
function already_found(df, v, prime; pmap_prime_tol = PMAP_PRIME_TOL)
    v === nothing && return false
    for o in eachrow(df)
        o.prime == prime || continue
        any(i -> norm([o.sec_y[i], o.sec_py[i]] .- v) < pmap_prime_tol, eachindex(o.sec_y)) &&
            return true
    end
    return false
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



"""
    Takes a bunch of seed and finds the periodic orbit to each seed. If the orbit is found, it is added to the df. The df is returned at the end.
"""
function sweep!(df::AbstractDataFrame, seeds, n, prm; pmap_root_tol = PMAP_ROOT_TOL, verbose = false)
    @showprogress dt = 1 desc = "E=$(round(prm.E; digits = 4)) n=$n" for v0 in seeds
        orb = try
            analyse_seed(v0, n, prm; str="seed", pmap_root_tol = pmap_root_tol)
        catch e
            verbose && @warn "seed failed" seed = v0 exception = e
            continue
        end
        orb===nothing || push!(df, orb)
    end
    return df
end



 

"""
finding three base orbits at energy
    return
    (;orb=(; df, timing, seeds), fig=(f, ax))
"""
function get_obrits_ABC(;
    p            = (1.0, 1.0, 1.0),
    Emin         = 0.01,
    nfast        = 1,              # crossings the dense integrator may take
    ndense       = 2,
    tmax         = 100_000.0,
    seeds = [[0.0,0.23], [0.2, 0.3], [0.3, 0.0]],
    orbit_str = ["A", "B", "C"],
    cmap      = COLOR_SCHEME,
    display_figure=false,
    show_figure=true
    )
    #seeds =  section_grid(E, p; ny = 3, npy = 3, margin = 0.03)
    any(!,[in_section(vi, Emin, p) for vi in seeds]) && error("some seeds are outside the energy boundary")

    
    prm = SectionParams(Emin, p, :both; tmax, nfast, ndense, save_everystep = false, save_start = false)

  

    df  = orbit_table()
    timing = DataFrame(E = Float64[], n = Int[], nseeds = Int[],
                       found = Int[], secs = Float64[])

    before = nrow(df)
    secs = @elapsed for i in axes(seeds, 1)
        v0 = seeds[i]

        sol = analyse_seed(v0, nfast, prm;str=orbit_str[i])
        sol !== nothing && push!(df, sol)
    end 
    push!(timing, (Emin, nfast, length(seeds), nrow(df) - before, secs))

    if display_figure  # surface section of the 
        cABC = [pick_color(i,cmap) for i in 1:3]
        der_rand =  boundary(Emin, p)
        f = Figure(size=(1400,900))
        ax = Axis(f[1, 1], xlabel = "y", ylabel = "py", title = "E=$(round(Emin, digits=4)) orbtis A, B, C")
        vs = Point2f.([o.v[1] for o in eachrow(df)], [o.v[2] for o in eachrow(df)])
        scatter!(ax, vs, color = cABC, markersize = 6, alpha = 1)
        annotation = ["$(orbit_str[i]), T=$(round(o.T, digits=3))" for (i, o) in enumerate(eachrow(df))]
        text!(ax, vs, text=annotation, fontsize = 14, color = INK, align = (:left, :bottom), offset = (5, 5))

        scatter!(ax, Point2f.(vs), markersize = 5, color = INK, alpha = 0.7)
        scatter!(ax, der_rand, markersize = 3, color = INK)
        show_figure && display(f)
        return (;orb=(; df, timing, seeds), fig=(f, ax))
    else
        return (;orb=(; df, timing, seeds), fig=nothing)
    end
end




"""
    follow the orbits A, B, C from low energies upward. The goal is to trace the Monodromy matrix as a function of energy. 
    Start at Energy E=0.11
    return orbs0
"""
function follow_ABC!(orbs0, Es; 
                 p=(1.0,1.0,1.0), tmax = 100_000.0, nfast = 1, ndense = 1, verbose = false)
    prm = SectionParams(orbs0.E[1], p, :both; tmax, nfast, ndense, save_everystep = false, save_start = false)

    for o in collect(eachrow(orbs0))
        v0          = o.v
        orbit_label = o.str
        df          = orbit_table()
        misses      = 0
        @showprogress dt=1 desc="orbit $orbit_label" for e in Es
            set_energy!(prm, e)
            sol = try
                analyse_seed(v0, nfast, prm; str = orbit_label)
            catch err
                verbose && @warn "step failed" orbit=orbit_label E=e exception=err
                nothing
            end
            if sol === nothing
                misses += 1
                misses ≥ 3 && (verbose && @info "orbit $orbit_label lost near E=$e"; break)
                continue
            end
            misses = 0
            push!(df, sol)
            v0 = sol.v
        end
        append!(orbs0, df)
    end
    return (; orb = orbs0)
end

"""
    return (;M, check)

    This takes care of (NaN or Inf) ∈ M
"""
function monodrome(orbits; p = (1.0,1.0,1.0), verbose=false)
    M    = Vector{Matrix{Float64}}(undef, nrow(orbits))
    check = falses(nrow(orbits))

    @showprogress dt=1 desc="monodromy" for (i, o) in enumerate(eachrow(orbits))
        Mat = try
            monodromy(lift(o.v, o.E, p), o.T; p=p)
        catch err
            verbose && @warn "monodromy failed" E=o.E str=o.str exception=err
            fill(NaN, 4, 4)
        end
        sane     = !(any(isnan, Mat) || any(isinf, Mat)) && abs(det(Mat) - 1.0) < 0.01
        M[i]     = sane ? Mat : fill(NaN, 4, 4)
        check[i] = sane    
    end

    return (; M, check)
end


"
appends .M and .check to orbits

    orbits.mono.M:: Matrix
    orbits.mono.check:: Bool(sensible matrix)
"
function append_monodrome(orbits::DataFrame)
    orb = copy(orbits)
    mono = monodrome(orbits)
    orb.M     = mono.M
    orb.check = mono.check
    return orb
end


function eigen_tr(M)
    τ = tr(M)
    p= (τ-2)
    return p/2 .+ (1,-1) .* sqrt(Complex((p/2)^2-1))
end

"
    M -> Monodromy matrix
    λ -> max(eigenvalues(M))
    T -> period of periodic orbit
    returns |α|T:: lyapunof exponent to λ= exp(|α|T)
"
get_lyapunov(M::Matrix, T)   = get_lyapunov(eigvals(M), T)
get_lyapunov(λs::AbstractVector, T) = log(maximum(abs.(λs))) 

get_phase(M::Matrix)          = get_phase(eigen_tr(M))
get_phase(λs)                 =  abs(angle(λs[1]))


function ABC_energy_trace(;nup=5000,ndown=5000)
    p            = (1.0, 1.0, 1.0)
    E_fix        = 0.11  # this is the one I fixed
    E_max        = 1.0
    Es_up        = collect(range(E_fix, E_max, nup))[2:end]
    E_min        = 0.01
    Es_down      = sort(collect(range(E_min, E_fix, ndown))[1:end-1]; rev=true)
    nfast        = 1              # crossings the dense integrator may take
    ndense       = 2
    tmax         = 100_000.0
    seeds = [[0.0,0.23], [0.2, 0.3], [0.3, 0.0]]
    orbit_str = ["A", "B", "C"]

    # starting point to follow A, B and C orbit
    res = get_obrits_ABC(p = p, 
                        Emin = E_fix,
                        nfast = nfast, 
                        ndense = ndense, 
                        tmax = tmax,
                        seeds = seeds,
                        orbit_str = orbit_str, 
                        display_figure = false)

    orb_ABC_up, orb_ABC_down=copy(res.orb.df), copy(res.orb.df)  # copy for clean split of variables 


    follow_ABC!(orb_ABC_down,   Es_down)  # follow to lower energy
    follow_ABC!(orb_ABC_up,     Es_up)    # follow to higher energy

    all_ABC = vcat(orb_ABC_down,orb_ABC_up)
    sort!(all_ABC,[:E, :T])

    return (; all_ABC)
end

"""
    return (;v0, E_bif, T_bif)

    lable = i.e. "A" / "B" / "C" etc...
    df::DataFrame that contains orbits
    order = Θ=p/order * 2π

Where E_bif is the Energy where the bifurcation takes place and T_bif that the orbit  has when bifurcation happens
"""
function get_bifur_E_T(label, df, order)
    # per orbit
    sub        = filter(o -> o.str == label, df)
    Ms, check  = ("M" ∈ names(sub) && "check" ∈ names(sub)) ? (sub.M, sub.check) : monodrome(sub)
    v0s        = sub.v[check]
    Es         = sub.E[check]
    Ts         = sub.T[check] 
    eigs       = [eigen_tr(M) for M in Ms[check]]
    θs         = [get_phase(collect(λs)) for λs in eigs] ./2pi
    is_bifurc = θs .> 1/order   # obtain point of bifurcation 0->1 or 1->0
    bifurc_exists = any(is_bifurc)
    # obtain period of bifurcation orbit
    v0, E_bif, T_bif = bifurc_exists ? (first(v0s[is_bifurc]), first(Es[is_bifurc]), order*first(Ts[is_bifurc])) : (fill(NaN, 2), NaN, NaN)
    isnan(E_bif) && println("="^72, "\nFor orbit $label, no bifurcation of order $order found!\n", "="^72)
    return (;v0, E_bif, T_bif)
end




"""
    return orbits' (that bifurcates with given 'order' of bifurcation. )

Finds bifurcation at Energy E(Θ) such that Θ ≈ (n/m *2 π), m=order.  
"""
function get_bifur_orbit(orbits::DataFrame, order::Int; labels=["A", "B", "C"],p =(1.0,1.0,1.0))
    df = sort(copy(orbits), [:E])
    bi_orbs = orbit_table()
    # get monodromy matrix of orbits
    orbA = get_bifur_E_T(labels[1], df, order)
    isnan(orbA.E_bif) && error("no order-$order bifurcation found on branch $(labels[1])")
    # get ready to find the orbit that intersects 'order'-times with sirface of section
    prm = SectionParams(orbA.E_bif, p, :both; tmax=(order + 1) * orbA.T_bif, nfast=order, ndense=order+1, save_everystep = false, save_start = false)
    if true
        y0, py0 = orbA.v0
        npoint  = 200
        yzoom   = 0.001
        pyzoom  = 0.02                       # box size around the parent periodic orbit

        y_grid  = range(y0 - yzoom,  y0 + yzoom,  npoint)
        py_grid = range(py0 - pyzoom, py0 + pyzoom, npoint)

        set_energy!(prm, orbA.E_bif)  # past the bifurcation -- AT E_bif the satellite
                                            # has zero amplitude, nothing to see there
        Z = [norm(Fres_safe([y, py], order, prm)) for y in y_grid, py in py_grid]

        f  = Figure(size = (1300, 900))
        ax = Axis(f[1, 1], xlabel = L"y", ylabel = L"p_y",
                title = "log10|T^$order(v)-v| around orbit A, E=$(round(prm.E, digits=4))")
        contour!(ax, y_grid, py_grid, log10.(Z .+ 1e-300))   # floor -- log10(0) = -Inf otherwise
        scatter!(ax, y0, py0, color = INK, markersize = 8)
        display(f)
    end
    function find_bif_orbit(v0, n, prm; str = "seed",
                        pmap_root_tol = PMAP_ROOT_TOL, pmap_prime_tol= PMAP_PRIME_TOL, maxiters=300)
        res = solve_orbit(v0, n, prm, tol=pmap_root_tol, maxiters=maxiters)
        res.converged || return nothing

        min_period = minPeriodicity(res.v, prm; pmap_prime_tol = pmap_prime_tol, search = 40)
        min_period.Nperiod === nothing && return nothing
        T = min_period.Tperiod
        return (; E = prm.E, v = res.v, str, T), min_period
    end   
    orb_bi, min_period  = find_bif_orbit(orbA.v0, order, prm; str="A.$order")
    orb_bi === nothing ?  error("The orbit sought for is not found!") : push!(bi_orbs, orb_bi)
    println("comparing periods; A: T=$(orbA.T_bif) , A.4 T= $(orb_bi.T)")
    println(bi_orbs)
    Mbi = monodrome(bi_orbs)
    println(typeof(Mbi.M[1]))
    L_Mbi = get_lyapunov(Mbi.M[1], orb_bi.T)
    println("Lyapunof exponent $L_Mbi")
    if false 
        set_style!(:dark)
        der_rand =  boundary(orbA.E_bif, p)
        v0       = Point2f(orbA.v0)
        v_bi     = Point2f(orb_bi.v)
        vpMapbi  = Point2f.(min_period.pMap)
        f = Figure(size=(1300, 900))
        ax= Axis(f[1,1], ylabel=L"p_y", xlabel=L"y", title="Surface of section at Bifurcation Energy $(orbA.E_bif)")
        scatter!(ax, v0, color = pick_color(1), markersize = 6, alpha = 1)
        map(v0i_bi->scatter!(ax, v0i_bi, color = pick_color(3), markersize = 6, alpha = 1), vpMapbi)
        scatter!(ax, v_bi, color = pick_color(2), markersize = 6, alpha = 1)
    
        annotation = ["A",orb_bi.str, "period"]
        text!(ax, [v0, v_bi, first(vpMapbi)], text=annotation, fontsize = 14, color = INK, align = (:left, :bottom), offset = (5, 5))
        scatter!(ax, der_rand, markersize = 3, color = INK)
        display(f)
    end

    return orb_bi
end



function graphs(res::DataFrame; fsize=(1400,900), cmap=COLOR_SCHEME)
    df = sort(copy(res), [:E])
    println("="^72)


    function fig_for(label, df)
        sub        = filter(o -> o.str == label, df)

        Ms, check  = ("M" ∈ names(sub) && "check" ∈ names(sub)) ? (sub.M, sub.check) : monodrome(sub)
        Es         = sub.E[check]
        eigmat     = reduce(hcat, eigvals.(Ms[check]))'   # rows = energy, cols = branch j=1:4
        eigs       = [eigen_tr(M) for M in Ms[check]]
        lyp        = map((λs,T) -> get_lyapunov(collect(λs),T), eigs, sub.T[check])
        traces      = [tr(m) for m in Ms[check]]
        θs         = [get_phase(collect(λs)) for λs in eigs] ./2pi

        fλ  = Figure(size = fsize)        
        fΘ  = Figure(size = fsize)
        fLy = Figure(size = fsize)
        fTr = Figure(size = fsize)      

        # Lyapunov exponents of (monodrome) Map 
        axLy = Axis(fLy[1,1], xlabel = "E", ylabel = L"|α|T", title = "orbit $label")  
        # trace of Monodromy
        axTr = Axis(fTr[1,1], ylabel = "tr(M)")
        hlines!(axLy,[1.0], color=INK, linewidth=5, linestyle = :dash)
        lines!(axLy, Es, lyp)
        lines!(axTr, Es, traces)

        # Θ plot -> for theta = p//q, p,q∈N, bifurcation happens
        pdivq = [0,1/2, 1/3, 1/4, 1/5, 2/5, 1/6]
        yticks= ["0",L"\frac{1}{2}", L"\frac{1}{3}", L"\frac{1}{4}", L"\frac{1}{5}", L"\frac{2}{5}", L"\frac{1}{6}"]
        axΘ = Axis(fΘ[1,1], xlabel = "E", ylabel = "Θ/2π", title = "orbit $label", yticks=(pdivq,yticks), yticklabelcolor = INK, yaxisposition = :left)
        lines!(axΘ, Es, θs)
        hlines!(axΘ, pdivq[2:end], linestyle=:dot, color=INK)
        label == "A" && vlines!(axΘ, get_bifur_E_T(label,sub, 6).E_bif, color=INK, linestyle=:dot)
        
        # eigenvalues λ, 1/λ
        axλ = Axis(fλ[1,1], xlabel = "E", ylabel = "Re(λ) and Im(λ)", title = "orbit $label")


        lines!(axλ, Es, real.(first.(eigs)); color = pick_color(1, cmap), linewidth = 2, linestyle = :solid)
        lines!(axλ, Es, imag.(first.(eigs)); color = pick_color(1, cmap), linewidth = 2, linestyle = :dash)
        lines!(axλ, Es, real.(last.(eigs)); color = pick_color(2, cmap), linewidth = 2, linestyle = :solid)
        lines!(axλ, Es, imag.(last.(eigs)); color = pick_color(2, cmap), linewidth = 2, linestyle = :dash)
            # text!(ax, Es[end], real(eigmat[end, j]); text = "λ$j",
            #       color = pick_color(j, cmap), align = (:left, :center), offset = (6, 0), fontsize = 13)
 

        style_elems  = [LineElement(color = :gray70, linestyle = :solid, linewidth = 2),
                        LineElement(color = :gray70, linestyle = :dash,  linewidth = 2)]
        branch_elems = [LineElement(color = pick_color(j, cmap), linewidth = 2) for j in axes(eigmat, 2)]

        Legend(fλ[1,2],
               [style_elems, branch_elems],
               [["Re(λ)", "Im(λ)"], ["λ$j" for j in [1,2]]],
               ["Component", "Branch"])

        (; fλ, axλ, fLy, axLy, fΘ, axΘ, fTr, axTr)
    end

    fA, fB, fC = fig_for("A",df), fig_for("B",df), fig_for("C",df)
    return (; fA , fB , fC)
end




# res = ABC_energy_trace(nup=5000, ndown=5000).all_ABC
# orbits = append_monodrome(res)
# savename = "orbitsABC.jld2"

# JLD2.save_object(joinpath(SAVE_DATA_DIR, savename), orbits)
orbits = JLD2.load_object(joinpath(SAVE_DATA_DIR, savename))

# looking at the contour around orbit that bifurcates
p=(1.0,1.0,1.0)  #  parameters
order=4          #  bifurcation order


orbA = get_bifur_E_T("A", orbits, order)

isnan(orbA.E_bif) && error("no order-$order bifurcation found on branch $(labels[1])")
    # get ready to find the orbit that intersects 'order'-times with sirface of section
prm = SectionParams(orbA.E_bif, p, :both; tmax=(order + 1) * orbA.T_bif,nfast=order, ndense=order+1, save_everystep = false, save_start = false)
if true
    set_style!(:print)
    y0, py0 = orbA.v0
    npoint  = 1000
    yzoom   = 0.0075
    pyzoom  = 0.003                       # box size around the parent periodic orbit

    y_grid  = range(y0 - yzoom,  y0 + yzoom,  npoint)
    py_grid = range(py0 - pyzoom, py0 + pyzoom, npoint)

    set_energy!(prm, orbA.E_bif)  # past the bifurcation -- AT E_bif the satellite

    Z = [norm(Fres_safe([y, py], order, prm)) for y in y_grid, py in py_grid]

    f  = Figure(size = (1300, 900))
    ax = Axis(f[1, 1], xlabel = L"y", ylabel = L"p_y",
                title = "After bifurcation order $order @ E=$(round(prm.E, digits=4))")
    contour!(ax, y_grid, py_grid, log10.(Z .+ 1e-300))   # floor -- log10(0) = -Inf otherwise
    scatter!(ax, y0, py0, color = INK, markersize = 8)

    
    save(joinpath(FIG_DIR, "contour_at_bif_order$order.png"), f; px_per_unit = 2)
    display(f)
end





set_style!(:dark)  # :print
figs = graphs(orbits)


# display(figs.fA.fλ)
# display(figs.fB.fλ)
# display(figs.fC.fλ)

# display(figs.fA.fLy)
# display(figs.fB.fLy)
# display(figs.fC.fLy)

display(figs.fA.fΘ)
# display(figs.fB.fΘ)
# display(figs.fC.fΘ)


# xlims!(figs.fB.ax, -0.001, 0.1)
# ylims!(figs.fB.ax, -0.5,3)
# save(joinpath(FIG_DIR, "orbitA/Aeigenvalues-vs-E.png"), figs.fA.fig; px_per_unit = 2)
