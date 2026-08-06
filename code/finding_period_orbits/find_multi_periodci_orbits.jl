import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2
include("../styles/makie_theme.jl")
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
    reinit!(pmap, [1e-10, v[1], px, v[2]])
    trace = zeros(2, n)
    for i in 1:n
        step!(pmap, 1)
        u = current_state(pmap)
        #println(u)
        trace[:, i] = [u[2], u[4]]
    end
    return trace
end



Fres(pmap, v, n, E, p) = Tn(pmap, v, n, E, p) .- v     # residual, roots = periodic orbits

# --- Newton iteration ---

function find_orbit(pmap, v0, n, E, p;
                    N_max = 50, d = 1e-6, damping = 1.0, tol = 1e-9)
    v       = collect(float.(v0))
    history = zeros(2, N_max)
    in_section(v, E, p) ||
        return (v = v, converged = false, history = history[:, 1:0], resnorm = Inf)

    JT = zeros(2, 2)
    rn = Inf
    for i in 1:N_max
        history[:, i] = v
        r  = Fres(pmap, v, n, E, p)
        rn = norm(r)
        rn < tol && return (v = v, converged = true, history = history[:, 1:i], resnorm = rn, comment="$rn < $tol")

        JT[:, 1] = (Fres(pmap, [v[1]+d, v[2]], n, E, p) .- Fres(pmap, [v[1]-d, v[2]], n, E, p)) ./ (2d)
        JT[:, 2] = (Fres(pmap, [v[1], v[2]+d], n, E, p) .- Fres(pmap, [v[1], v[2]-d], n, E, p)) ./ (2d)
        for (i, f) in enumerate(JT)
            if abs(f) < 1e-11
                println("derivative is zero at $i")
            end 
        end
        vnew = v .- damping .* (JT \ r)
        in_section(vnew, E, p) ||
            return (v = v, converged = false, history = history[:, 1:i], resnorm = rn, comment="not in section")
        v = vnew
    end
    (v = v, converged = false, history = history, resnorm = rn, comment="$N_max iteratrions passed")
end

function init_grid(E, dy0, dpy0, y_max, py_max, p)
    y_min, py_min
    y0 = collect(y_min[1]:dy0:y_max[2])
    py0 = collect(-py_min:dpy0:py_max)
    v = zeros(Float64, 2)
    for y0_i, i in enumerate(y0)
        for py0_j, j in enumerate(py0)
            v = [y0_i, py0]
            if in_section(v, E, p)
                y0[i], py0[i]= y0_i, py0_i
            else
                y0[i], py0[i]= NaN, NaN
            end
        end
    end
    print("lengths of\\y0: $(length(y0))\\py0: $(length(py0))\\ \\$(length(y0)*length(py0))")
    y0, py0
end

# --- run ---
n     = [1,2,3,4,5,6,7,8,9]   # period of orbit sought
dy0    = 0.05
dpy0   = 0.05
y0, py0 = init_grid(E0, dy0, dpy0, param)

for n_i, l in enumerate(n)
    for y0_i, i in enumerate(y0)
        for py0_j, j in enumerate(py0)
            if y0_i!=NaN && py0_j!=NaN
                u0    = [x0, y0_i, px_from_E(y0_i, py0, E0, param), py0_j]
                ds   = CoupledODEs(equations!, u0, param;
                                diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9))
                pmap = PoincareMap(ds, plane; direction = -1,
                                rootkw = (xrtol = 1e-12, xatol = 1e-12))
                res = find_orbit(pmap, [y0, py0], n, E0, param, N_max=1000, damping=0.1,tol=1e-9)
                println(res.converged ? "converged: v = $(res.v), distance $(res.resnorm), comment: $(res.comment)" : "no convergence, $(res.resnorm), comment: $(res.comment)")
                
                name = "period_of$(n)-aty$y0-py$py0.png"
                outfig = joinpath(FIG_DIR, name)

                with_theme(QUARTO_THEME) do
                    fig = Figure(size = (1920, 1200))
                    ax  = Axis(fig[1, 1], xlabel = "y", ylabel = "p_y")
                    scatter!(ax, y_all, py_all, colormap=:viridis, markersize=1.5)
                    scatter!(ax, boundary, color=C_CREAM)
                    scatterlines!(ax, res.history[1, :], color=C_PURPLE, res.history[2, :])
                    scatter!(ax, [res.v[1]], [res.v[2]], marker=:cross, markersize = 20, color=C_CREAM)
                    save_screen = GLMakie.Screen(; visible=false)
                    scatter!(ax, y_sec, py_sec, color=C_TEAL)
                    display(save_screen, fig)
                    save(outfig, fig, px_per_unit=1)
                    close(save_screen)
                    #display(GLMakie.Screen(), fig)
                end  # quarto theme
            else
                nothing
            end
        end
    end
end

u0     = [0.0, res.v[1], px_from_E(res.v[1], res.v[2], E0, param), res.v[2]]
y_sec, py_sec = Float64[], Float64[]
cb = HenonHeiles.section_callback(y_sec, py_sec)
s = HenonHeiles.solve_trajectory(u0, param, (0.0, T), dt; callback=cb)



# --- plot ---
