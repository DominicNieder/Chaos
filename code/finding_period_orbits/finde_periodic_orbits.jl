import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2,
      NonlinearSolve, ADTypes
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

# --- background section ---
data = load(data_file, "results")
y_all, py_all = Float32[], Float32[]
for d in data
    append!(y_all,  d.sec_y)
    append!(py_all, d.sec_py)
end

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


px2(y, py, E, p)    = 2p[2] * (E - HenonHeiles.potential(0.0, y, p)) - py^2
in_section(v, E, p) = px2(v[1], v[2], E, p) > 0




"Lift (y, py) to a 4D state. Offset follows sign(px) -> no phantom t=0 crossing."
function lift(v, E, p; sgn = +1)
    a = px2(v[1], v[2], E, p)
    a <= 0 && return nothing
    return [EPS_OFF, v[1], sgn * sqrt(a), v[2]]
end

function section_trj(v, prm)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")

    y_sec, py_sec, t_sec = Float64[], Float64[], Float64[]
    cb = HenonHeiles.section_callback(y_sec, py_sec, section_t=t_sec)
    prob = ODEProblem(HenonHeiles.equations!, u0, (0.0, 2000.0), prm.p)
    sol = solve(prob, Vern9(), 
        abstol=1e-14, 
        reltol=1e-14, 
        saveat=dt, 
        callback=cb
        )

    return  sol, permutedims([y_sec py_sec]), t_sec
end

"returns named tuple (;traj, pMap, Nperiod, Tperiod)"
function minPeriodicity(v, prm; tol = 1e-8)
    sol, trace, ts = section_trj(v, prm)
    for i in axes(trace, 2)
        if (norm(v .- trace[:, i]) < tol)   
            k = searchsortedfirst(sol.t, ts[i])
            return (;traj=sol.u[1:k], pMap=trace[:,1:i], Nperiod=i, Tperiod=ts[i])
        end
    end
    return (;traj=sol.u, pMap=trace, Nperiod=nothing, Tperiod=nothing)
end

function creat_integrator(p; tmax = 2000.0, nmax = Ref(1))
    y, py, ts = Float64[], Float64[], Float64[]
    condition(u, t, integ) = u[1]
    affect!(integ) = begin
        push!(y, integ.u[2]); push!(py, integ.u[4]); push!(ts, integ.t)
        length(y) ≥ nmax[] && terminate!(integ)
    end
    cb = ContinuousCallback(condition, affect!, nothing; abstol = 1e-13)
    prob = ODEProblem(HenonHeiles.equations!, zeros(4), (0.0, tmax), p)
    integ = init(prob, Vern9(); abstol = 1e-14, reltol = 1e-14,
                 save_everystep = false, save_start = false, callback = cb)
    return (; integ, y, py, ts, nmax)
end

"prm= (; integ, n, E, p)"
function T(v,n, prm)
    u0 = lift(v, prm.E, prm.p)

    u0 === nothing && error("point $v outside energy boundary")
    # init problem 
    empty!(prm.y); empty!(prm.py); empty!(prm.ts)
    reinit!(prm.integ, u0)
    solve!(prm.integ)

    length(prm.y) < n &&
        error("only $(length(prm.y)) crossings in t < $(prm.tmax) (need $(n))")

    return [prm.y[n], prm.py[n]]
end


Fres(v, n, prm) = T(v, n, prm) - v

function Fres_safe(v, n, prm)
    a = px2(v[1], v[2], prm.E, prm.p)
    a <= 0 && return fill(1.0 + 100*sqrt(-a), 2)   # penalty scales with violation
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
            error("both probes outside boundary at $v — step too large or v on the edge")
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

function find_orbit(v0, prm;
                    N_max = 100, d = 1e-7, tol = 1e-11, 
                    max_backtrack = 30, dmax = 0.05, verbose = false)


    v    = collect(float.(v0))
    hist = zeros(2, N_max)
    J    = zeros(2, 2)

    fail(i, rn, msg) = (v = v, DT = nothing, eigen = nothing, converged = false,
                        history = hist[:, 1:max(i, 0)], resnorm = rn, comment = msg)

    in_section(v, prm.E, prm.p) || return fail(0, Inf, "seed outside energy boundary")

    r  = Fres(v, prm.n, prm)
    rn = norm(r)

    for i in 1:N_max
        hist[:, i] = v

        if rn < tol
            jacobian!(J, v, prm.n, prm, d)          # evaluated AT the root
            DT = J + I
            return (v = v, DT = DT, eigen = eigvals(DT), converged = true,
                    history = hist[:, 1:i], resnorm = rn,
                    comment = "|r| = $rn  det(DT) = $(det(DT))")
        end

        jacobian!(J, v, prm.n, prm, d)
        κ = cond(J)
        κ > 1e12 && @warn "ill-conditioned Jacobian" iterate=i cond=κ

        step = J \ r
        sn   = norm(step)
        sn > dmax && (step .*= dmax / sn)     # cap: keeps Newton local

        λ, accepted = 1.0, false              # MUST start false
        for _ in 1:max_backtrack
            vnew = v .- λ .* step
            rnew = in_section(vnew, prm.E, prm.p) ? Fres(vnew, prm.n, prm) : nothing
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

function solve_orbit(v0, prm; tol = 1e-11, maxiters = 300)
    v = collect(float.(v0))
    in_section(v, prm.E, prm.p) ||
        return (v = v, DT = nothing, eigen = nothing, converged = false,
                resnorm = Inf, history = zeros(2,0), comment = "seed outside boundary")

    prob = NonlinearProblem((w, q) -> Fres_safe(w, q.n, q), v, prm)    
    sol  = solve(prob, TrustRegion(autodiff = AutoFiniteDiff(fdtype = Val(:central)));
                 abstol = tol, maxiters)

    ok = SciMLBase.successful_retcode(sol) && in_section(sol.u, prm.E, prm.p)
    ok || return (v = sol.u, DT = nothing, eigen = nothing, converged = false,
                  resnorm = norm(sol.resid), history = zeros(2,0),
                  comment = "$(sol.retcode)")

    DT = jacobian!(zeros(2,2), sol.u, prm.n, prm, 1e-7) + I
    return (v = sol.u, DT = DT, eigen = eigvals(DT), converged = true,
            resnorm = norm(sol.resid), history = zeros(2,0),
            comment = "$(sol.retcode)  det(DT) = $(det(DT))")
end

"Is the orbit invariant under S(y,py) = (y,-py)?"
function sym_py(sec; tol = 1e-6)
    Ssec = sec .* [1.0, -1.0]
    for j in axes(Ssec, 2)
        any(i -> norm(Ssec[:, j] - sec[:, i]) < tol, axes(sec, 2)) || return false
    end
    return true
end

function Orbit_finder(v0, n, E;
    psection_bg = (y_all, py_all), tmax=20_000.0)
    println("=== NEW SIM ===")


    p       =   (1.0, 1.0, 1.0)

    bff     =   creat_integrator(p; tmax=tmax, nmax=Ref(n))

    prm = (; integ = bff.integ, y = bff.y, py = bff.py, ts = bff.ts,
                 E, p, n, tmax)

    res     =   find_orbit(v0, prm)
    orbit   =   minPeriodicity(res.v, prm; tol = 1e-8)
    u, sec, prime, t_sec    =  (orbit.traj, orbit.pMap, orbit.Nperiod, orbit.Tperiod)   # (;traj=sol.u, pMap=trace, Nperiod=nothing, Tperiod=nothing)
    

    # res.DT =DT
    if res.converged 
        DT=get_DT(res.v, prime, prm)

        println("found periodic orbit... symmetrizing...")
        if !sym_py(sec) # this returns the symmetrized section map
            sym_v = copy(res.v) .* [1,-1]
            sym_orbit= minPeriodicity(sym_v, prm; tol=1e-8)
            sym_u, sym_sec, sym_prime, sym_t_sec = (sym_orbit.traj, sym_orbit.pMap, sym_orbit.Nperiod, sym_orbit.Tperiod)
            DT_sym=get_DT(sym_v, prime, prm)
            println("symmertrized orbit!")
            println("comparing DT_sym:\n$(DT_sym)")
            println("DT:\n$(res.DT)")
            println("sym prime: $sym_prime")
        else
            sym_orbit = orbit
            sym_u, sym_sec = u, sec
            DT_sym = DT
            println("orbit is S-invariant: S(sec) = sec")
        end

        τ = tr(DT)
        t_sym=tr(DT_sym)
        kind = kind_index(τ; ε = 1e-6)         
        kind_sym=kind_index(t_sym; ε = 1e-6)
        println("CONVERGED")   
        println("   v         = $(res.v)")
        println("  |r|        = $(res.resnorm)")
        println("  det(DT)    = $(det(DT)), ($(det(DT_sym)))     (must be ~1)")
        println("  eigen      = $(res.eigen), ($(eigvals(DT_sym)))")
        println("  type       = $(KIND_LABEL[kind])   ($(KIND_LABEL[kind_sym]))")
        println("iterations n = $prime")
        println("orbit period = $t_sec")
    else
        println("did not converge: $(orbit.comment)")
    end
    


    println("=== plotting ===")

    
    x = [x[1] for x in u]
    y = [x[2] for x in u]

    # trajectory config space
    r      = range(-1.0, 1.0, length=120)
    levels = logrange(5.0*0.009, 6.9*0.089, 7)
    epot   = [HenonHeiles.potential(x,y, param) for x in r, y in r]

    fig_conf = Figure(size = (1400, 900))
    ax_conf  = Axis(fig_conf[1, 1], xlabel = "x", ylabel = "y",
                title = "config space, n = $prime", aspect = DataAspect())
    contour!(ax_conf, r, r, epot, labels=true, levels=levels, colormap=:hsv, colorscale=identity)
    lines!(ax_conf, x, y, color = C_PURPLE, label = "orbit, T=$t_sec")
    lines!(ax_conf, [x[1] for x in sym_u], [x[2] for x in sym_u], color = C_PURPLE, label = "sym. orbit, T=$t_sec")

    display(GLMakie.Screen(), fig_conf)

    fig_p = Figure(size=(1400,900))
    ax_p = Axis(fig_p[1, 1], xlabel = "y", ylabel = "p_y",
               title = "section, n = $n, prime = $(prime)")
    
    ybg, pybg = psection_bg
    y_max, py_max = HenonHeiles.section_boundary_ranges(E0, param, 120)
    boundary      = HenonHeiles.section_boundary(y_max, py_max) 
    scatter!(ax_p, ybg, pybg, color = (:salmon, 0.5), markersize = 1.5)
    scatter!(ax_p, boundary, color = C_CREAM, markersize = 4)
    scatterlines!(ax_p, res.history[1, :], res.history[2, :],
                    color = C_GOLD, markersize = 6, label = "Newton path")
    scatter!(ax_p, sec[1,:], sec[2,:], color= C_GREEN ,markersize = 5, label="periodic orbit")  # 
    scatter!(ax_p, sym_sec[1,:], sym_sec[2,:], color= C_TEAL ,markersize = 5, label="sym. periodic orbit") 
    display(GLMakie.Screen(), fig_p)
    return (;pfig=(fig_p,ax_p), cfig=(fig_conf,ax_conf), res=res, orbit=orbit, sym_orbit, prm=prm)
end


E0 = 0.1127
n  = 6
y0 = -0.012
py0= 0.05
v0 = [y0,py0]
sol=Orbit_finder(v0, n, E0)