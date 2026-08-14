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
E0 = 0.1127


px2(y, py, E, p)    = 2p[2] * (E - HenonHeiles.potential(0.0, y, p)) - py^2
in_section(v, E, p) = px2(v[1], v[2], E, p) > 0

"Lift (y, py) to a 4D state. Offset follows sign(px) -> no phantom t=0 crossing."
function lift(v, E, p; sgn = +1)
    a = px2(v[1], v[2], E, p)
    a <= 0 && return nothing
    return [sgn*EPS_OFF, v[1], sgn * sqrt(a), v[2]]
end

function section_trj(v, prm)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")

    y_sec, py_sec, t_sec = Float64[], Float64[], Float64[]
    cb = HenonHeiles.section_callback(y_sec, py_sec, section_t=t_sec)
    prob = ODEProblem(HenonHeiles.equations!, u0, (0, 2_000), prm.p)
    sol = solve(prob, Vern9(), 
        abstol=1e-14, 
        reltol=1e-14, 
        saveat=dt, 
        callback=cb
        )

    return sol.u, permutedims([y_sec py_sec]), t_sec
end

function minPeriodicity(v, prm; tol = 1e-8)
    u, trace, ts = section_trj(v, prm)
    for i in axes(trace, 2)
        norm(v .- trace[:, i]) < tol && return (u, trace, i, ts[i])
    end
    return (u, trace, nothing, nothing)
end


"prm= (; pmap, n, E, p, sgn)"
function Tn(v, prm)
    u0 = lift(v, prm.E, prm.p)
    tmax =2_000
    u0 === nothing && error("point $v outside energy boundary")
    # init problem 
    y_sec, py_sec = Float64[], Float64[]
    cb = HenonHeiles.section_callback(y_sec, py_sec)
    prob = ODEProblem(HenonHeiles.equations!, u0, (0, tmax), prm.p)
    sol = solve(prob, Vern9(), abstol=1e-9, 
        reltol=1e-9, 
        saveat=dt, 
        callback=cb
        )
    length(y_sec) < prm.n &&
        error("only $(length(y_sec)) crossings in t < $tmax (need $(prm.n))")

    return [y_sec[prm.n], py_sec[prm.n]]
end


Fres(v, prm) = Tn(v, prm) - v


function jacobian!(J, v, prm, d)
    E, p = prm.E, prm.p
    for (j, e) in enumerate(([d, 0.0], [0.0, d]))
        vp, vm = v .+ e, v .- e
        okp, okm = in_section(vp, E, p), in_section(vm, E, p)
        if okp && okm
            J[:, j] = (Fres(vp, prm) .- Fres(vm, prm)) ./ (2d)
        elseif okp
            J[:, j] = (Fres(vp, prm) .- Fres(v, prm)) ./ d       # forward
        elseif okm
            J[:, j] = (Fres(v,  prm) .- Fres(vm, prm)) ./ d      # backward
        else
            error("both probes outside boundary at $v — step too large or v on the edge")
        end
    end
    return J
end


function find_orbit(v0, n, E, p;
                    N_max = 100, d = 1e-7, tol = 1e-11, sgn = +1,
                    max_backtrack = 30, dmax = 0.05, verbose = false)

    prm  = (;n, E, p)
    v    = collect(float.(v0))
    hist = zeros(2, N_max)
    J    = zeros(2, 2)

    fail(i, rn, msg) = (v = v, DT = nothing, eigen = nothing, converged = false,
                        history = hist[:, 1:max(i, 0)], resnorm = rn, comment = msg)

    in_section(v, E, p) || return fail(0, Inf, "seed outside energy boundary")

    r  = Fres(v, prm)
    rn = norm(r)

    for i in 1:N_max
        hist[:, i] = v

        if rn < tol
            jacobian!(J, v, prm, d)          # evaluated AT the root
            DT = J + I
            return (v = v, DT = DT, eigen = eigvals(DT), converged = true,
                    history = hist[:, 1:i], resnorm = rn,
                    comment = "|r| = $rn  det(DT) = $(det(DT))")
        end

        jacobian!(J, v, prm, d)
        κ = cond(J)
        κ > 1e12 && @warn "ill-conditioned Jacobian" iterate=i cond=κ

        step = J \ r
        sn   = norm(step)
        sn > dmax && (step .*= dmax / sn)     # cap: keeps Newton local

        λ, accepted = 1.0, false              # MUST start false
        for _ in 1:max_backtrack
            vnew = v .- λ .* step
            rnew = in_section(vnew, E, p) ? Fres(vnew, prm) : nothing
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


function Orbit_finder(y0, py0, n, E;
    psection_bg = (y_all, py_all))


    v0= [y0, py0]
    p= (1.0, 1.0, 1.0)
    prm=(E=E, p=p, n=n)
    res = find_orbit(v0, n, E, p)

    if res.converged 
        println("CONVERGED   v = $(res.v)")
        println("  |r|        = $(res.resnorm)")
        println("  det(DT)    = $(det(res.DT))     (must be ~1)")
        println("  tr(DT)     = $(tr(res.DT))")
        println("  eigen      = $(res.eigen)")
    else
        println("did not converge: $(orbit.comment)")
    end
        # find period time and orbit and section trajectory
    (u, sec, prime, t_sec)=minPeriodicity(res.v, prm; tol = 1e-8)

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
    scatter!(ax_p, sec[:,1], sec[:,2], color= C_GREEN ,markersize = 5, label="periodic orbit")
    display(GLMakie.Screen(), fig_p)
    
end


n=7
y0= -0.075
py0=0.13
sol=Orbit_finder(y0, py0, n, E0)