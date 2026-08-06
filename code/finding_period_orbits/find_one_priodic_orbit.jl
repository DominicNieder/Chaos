import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2, NonlinearSolve
include("../styles/makie_theme.jl")
set_theme!(QUARTO_THEME)
include("../models/henon_heiles.jl")
using .HenonHeiles

CONFIG_DIR      = joinpath(@__DIR__, "../sim_config/henon_heiles.json")

DATA_DIR        = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
data_file       = joinpath(DATA_DIR, "E0.1127-T10000.0-py0.0-n256.jld2")

FIG_DIR         = joinpath(@__DIR__, "../../figures/henon-heiles/periodic-orbits/")



FIG_ORIENTATION  = joinpath(FIG_DIR, "orientation.json")
DATA_ORIENTATION = joinpath(DATA_DIR, "orientation.json")

configurations      = JSON3.read(read(CONFIG_DIR, String))
cfg                 = configurations.explore
param               = [Float64(cfg.a.value), Float64(cfg.m.value), Float64(cfg.w.value)]

T                   = 10^5
dt                  = Float64(cfg.dt.value)
x0                  = Float64(cfg.x0.value)


# --- phase map ---
data = load(data_file, "results")
y_all, py_all, c_all = Float32[], Float32[], Float32[]
for (i, d) in enumerate(data)
    append!(y_all, d.sec_y)
    append!(py_all, d.sec_py)
    append!(c_all, fill(Float32(i), length(d.sec_y)))
end

# --- plot Newton path over the section boundary ---
E0    = 0.1127
y_max, py_max = HenonHeiles.section_boundary_ranges(E0, param, 120)
boundary      = HenonHeiles.section_boundary(y_max, py_max)


# --- section map machinery ---

plane = (1, 0.0)                                      # x = 0

function px_from_E(y, py, E, p)
    arg = 2p[2]*(E - HenonHeiles.potential(0.0, y, p)) - py^2
    arg < 0 ? NaN : sqrt(arg)                          # NaN outside boundary
end

in_section(v, E, p) = !isnan(px_from_E(v[1], v[2], E, p))

"T^n final point: v = [y, py] -> [y', py'] after n crossings."
function Tn(pmap, v, n, E, p)
    px = px_from_E(v[1], v[2], E, p)
    isnan(px) && error("point $v outside energy boundary")
    reinit!(pmap, [1e-10, v[1], px, v[2]])
    step!(pmap, n)
    u = current_state(pmap)
    return [u[2], u[4]]
end
 
"All n crossings as a 2 x n matrix (for scans / period checking)."
function iterate_section(pmap, v, n, E, p)
    px = px_from_E(v[1], v[2], E, p)
    isnan(px) && error("point $v outside energy boundary")
    reinit!(pmap, [1e-15, v[1], px, v[2]])
    trace = zeros(2, n)
    for i in 1:n
        step!(pmap, 1)
        u = current_state(pmap)
        #println(u)
        trace[:, i] = [u[2], u[4]]
    end
    return trace
end


# Find residuals -> periodic orbits wth param= (;pmap, n, E, p)
Fres(v, p) = Tn(p.pmap, v, p.n, p.E, p.p) .- v     


function jacobian_get!(JT, pmap, v, n, E, p, d)
    JT[:, 1] = (Fres(pmap, [v[1]+d, v[2]], n, E, p) .- Fres(pmap, [v[1]-d, v[2]], n, E, p)) ./ (2d)
    JT[:, 2] = (Fres(pmap, [v[1], v[2]+d], n, E, p) .- Fres(pmap, [v[1], v[2]-d], n, E, p)) ./ (2d)
end
# --- Newton iteration method ---

get_eigenvalues(DT)= tr(DT)^2/4+sqrt(tr(DT)^2/4-det(DT)), tr(DT)^2/4-sqrt(tr(DT)^2/4-det(DT))

function find_orbit(pmap, v0, n, E, p;
                    N_max = 50, d = 1e-6, damping = 1.0, tol = 1e-9)
    v       = collect(float.(v0))
    history = zeros(2, N_max)
    in_section(v, E, p) ||
        return (v = v, DT = nothing, eigen=nothing converged = false, orbitTime=Inf history = history[:, 1:0], resnorm = Inf, comment="inital condition out of bounds!")

    JT = zeros(2, 2)
    rn = Inf
    for i in 1:N_max
        history[:, i] = v  # save section crossings to find roots
        r  = Fres(pmap, v, n, E, p)     # is sought to be zero
        rn = norm(r)                    # under this tolerance

        rn < tol && return (v = v, DT = JT+I, eigen=get_eigenvalues(JT+I), converged = true, orbitTime=current_time(pmap), history = history[:, 1:i], resnorm = rn, comment="$rn < $tol")

        jacobian_get!(JT, pmap, v, n, E, p, d)  #  jacobian to move toward Fres()=0

        for (indx_f, f) in enumerate(JT)
            if abs(f) < 1e-11
                println("derivative is zero at $indx_f")
            end 
        end

        vnew = v .- damping .* (JT \ r)  # iteration
        in_section(vnew, E, p) ||
            return (v = v, DT = nothing, eigen=nothing converged = false, orbitTime=Inf, history = history[:, 1:i], resnorm = rn, comment="not in section")  # check for bounds pf section
        v = vnew  # iterate
    end
    (v = v, DT = JT+ I, eigen=nothing , converged = false, orbitTime=Inf, history = history, resnorm = rn, comment="$N_max iteratrions passed")
end


function opt_perObits(pmap, u0, n, E, p;
    Nmax=200, tol=1e-12)
    v       = collect(float.(u0))
    param = (pmap=pmap, n=n, E=E, p=p)
    in_section(v, E, p) ||
        return (v = v, DT = nothing, eigen=nothing converged = false, orbitTime=Inf history = history[:, 1:0], resnorm = Inf, comment="inital condition out of bounds!")
    prob = NonlinearProblem(Fres, v, param)
    sol = solve(prob,
            TrustRegion(autodiff = AutoFiniteDiff(fdtype = Val(:central)));
            abstol   = tol,
            maxiters = Nmax)
end
# --- run ---
n     = 1   # period of orbit sought
y0    = -0.025
py0   = 0.128
u0    = [x0, y0, px_from_E(y0, py0, E0, param), py0]

ds   = CoupledODEs(equations!, u0, param;
                   diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9))
# I want to compare diratins
pmap = PoincareMap(ds, plane; direction = -1,
                   rootkw = (xrtol = 1e-12, xatol = 1e-12))
res = find_orbit(pmap, [y0, py0], n, E0, param, N_max=2000, damping=0.1,tol=1e-12)

pmap2 = PoincareMap(ds, plane; direction = +1,
                   rootkw = (xrtol = 1e-12, xatol = 1e-12))
res2 = find_orbit(pmap2, [y0, py0], n, E0, param, N_max=2000, damping=0.1,tol=1e-12)



println(res.converged ? "converged: v = $(res.v),time for reoccurence $(res.orbitTime). distance $(res.resnorm), comment: $(res.comment)" : "no convergence, $(res.resnorm), comment: $(res.comment)")
println(res2.converged ? "converged: v = $(res2.v), time for reoccurence $(res2.orbitTime), distance $(res2.resnorm), comment: $(res2.comment)" : "no convergence, $(res2.resnorm), comment: $(res2.comment)")



DT = res.DT
τ  = tr(DT)
@show det(DT), τ
println(abs(τ) < 2 ? "elliptic (stable)" :
        τ >  2      ? "hyperbolic (unstable)" :
                      "inverse hyperbolic (unstable, reflection)")
println("eigen values: $((τ/2)^2+sqrt((τ/2)^2-det(DT))), $((τ/2)^2-sqrt((τ/2)^2-det(DT)))")
# check for poincare secion that the found residual returns...
u01     = [1e-10, res.v[1], px_from_E(res.v[1], res.v[2], E0, param), res.v[2]]
y_sec, py_sec = Float64[], Float64[]
cb = HenonHeiles.section_callback(y_sec, py_sec)
s = HenonHeiles.solve_trajectory(u01, param, (0.0, 1000), dt; callback=cb)

u02     = [1e-10, res2.v[1], px_from_E(res2.v[1], res2.v[2], E0, param), res2.v[2]]
y_sec2, py_sec2 = Float64[], Float64[]
cb2 = HenonHeiles.section_callback(y_sec2, py_sec2)
s1 = HenonHeiles.solve_trajectory(u02, param, (0.0, 1000), dt; callback=cb2)

# contours of potential
r      = range(-1.0, 1.0, length=120)
levels = logrange(5.0*0.009, 6.9*0.089, 7)
epot   = [HenonHeiles.potential(x,y, param) for x in r, y in r]

fig2 = Figure(size = (1920, 1200))
ax2  = Axis(fig2[1, 1], xlabel = "X", ylabel = "y")
contour!(ax2, r, r, epot, labels=true, levels=levels, colormap=:hsv, colorscale=identity)
lines!(ax2, s1[1, :], s1[2, :], color = C_PURPLE)
lines!(ax2, s[1, :], s[2, :],color = C_TEAL)
display(fig2)


# --- plot ---
name = "period_of$(n)-aty$y0-py$py0.png"
outfig = joinpath(FIG_DIR, name)


fig = Figure(size = (1920, 1200))
ax  = Axis(fig[1, 1], xlabel = "y", ylabel = "p_y")
scatter!(ax, y_all, py_all, colormap=:viridis, markersize=1.5)
scatter!(ax, boundary, color=C_CREAM)
scatterlines!(ax, res.history[1, :], color=C_PURPLE, res.history[2, :])
scatterlines!(ax, res2.history[1, :], color=C_PURPLE, res2.history[2, :])

scatter!(ax, [res.v[1]], [res.v[2]], marker=:cross, markersize = 20, color=C_CREAM, label=-1)
scatter!(ax, [res2.v[1]], [res2.v[2]], marker=:cross, markersize = 20, color=:blue, label="+1")

save_screen = GLMakie.Screen(; visible=false)
scatter!(ax, y_sec, py_sec, color=C_TEAL)
scatter!(ax, y_sec2, py_sec2, color=C_GREEN)
display(save_screen, fig)
#save(outfig, fig, px_per_unit=1)
close(save_screen)
display(GLMakie.Screen(), fig)
