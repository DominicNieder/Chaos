using GLMakie

include("../models/henon_heiles.jl")
include("../analysis/section.jl")
include("../styles/makie_theme.jl")


with_theme(QUARTO_THEME) do
    # --- simulation variables ---
    E  = Observable(Float64(init_var[1]))
    x0 = Observable(Float64(init_var[2]))
    y0 = Observable(Float64(init_var[3]))
    py0= Observable(Float64(init_var[4]))
    T  = Observable(   Int64(num_int[1]))

    init = lift(E, x0, y0, py0, T) do e, x, y, py, T
        a, m, w = param
        v       = m*w^2/2*(x^2 + y^2) + a*(x^2*y - y^3/3)
        pmax    = sqrt(max(0.0, 2*e*m - v))
        py_c    = clamp(py, 0.0, pmax)
        diff    = 2e*m - m*w^2/2*(y^2 + x^2) - a*(x^2*y - y^3/3) - py_c^2
        px      = sqrt(max(0.0, diff))
        (x, y, px, py_c, T)
    end

    # --- derived trajectory ---
    sol = lift(init) do (x, y, px, py, T)
        solve_trajectory([x, y, px, py], tspan=(0.0, T), dt=num_int[2], p=param)
    end

    # --- (phase space) trajectory ---
    traj_x  = lift(s -> [u[1] for u in s.u], sol)
    traj_y  = lift(s -> [u[2] for u in s.u], sol)
    y_sec, py_sec = lift(s -> surface_section_x0(s), sol) |> s -> (lift(x->x[1],s), lift(x->x[2],s))

    # --- layout ---
    fig = Figure(size=(1920 , 1200))
    # real space trajectory
    ax1 = Axis(fig[1,1], title="Space (x, y)", xlabel="x", ylabel="y")
    xlims!(ax1, -1.0, 1.0)
    ylims!(ax1, -1.0, 1.0)
    lines!(ax1, traj_x, traj_y, linewidth=1.0)
    # plotting the contour of potential into fig 1
    r = range(-1.0, 1.0, length = 120)
    levels = logrange(5.0*0.009,6.9*0.089,7)
    epot = [potential(x,y; p=param) for x in r, y in r]
    contour!(ax1, r, r, epot, labels=true, levels=levels, colormap=:hsv, colorscale=identity)

    # phase space trajectory
    ax2 = Axis(fig[1,2], title="Section (x=0)", xlabel="y", ylabel="p_y")
    xlims!(ax2, -1.0, 1.0)
    ylims!(ax2, -1.0, 1.0)
    scatter!(ax2, y_sec, py_sec, markersize=5)

    # x(t)
    ax3 = Axis(fig[2,1:2], title="x(t)", xlabel="t", ylabel="x")
    lines!(ax3, lift(s->s.t, sol), traj_x, linewidth=0.5)
    ylims!(ax3,  -1.0, 1.0)
    on(T; update=true) do t
        xlims!(ax3, 0.0, t)
    end

    # --- controls ---
    sg = SliderGrid(fig[3, 1:2],
        (label="Energy E",  range=0.01:0.0001:0.3,   startvalue=0.083),
        (label="x₀",        range=-1.0:0.01:1.0,    startvalue=0.0),
        (label="y₀",        range=-0.8:0.01:1.2,    startvalue=0.0),
        (label="py₀",       range=0.0:0.01:1.0,     startvalue=0.0),
        (label="T",       range=100:50:10000,     startvalue=num_int[1])
    )
    on(sg.sliders[1].value) do v; E[]   = v; end
    on(sg.sliders[2].value) do v; x0[]  = v; end
    on(sg.sliders[3].value) do v; y0[]  = v; end
    on(sg.sliders[4].value) do v; py0[] = v; end
    on(sg.sliders[5].value) do v; T[]   = v; end


    Label(fig[4, 1:2], lift(v -> "effective py₀ = $(round(v[4], digits=4))", init))

    display(GLMakie.Screen(), fig)
end