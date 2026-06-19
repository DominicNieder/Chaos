using GLMakie

include("../models/henon_heiles.jl")
include("../analysis/section.jl")
include("../styles/makie_theme.jl")

# --- functions ---
function save_axis(ax, traj_observables...; fname)
    fig_out = Figure(size=(800, 700))
    ax_out  = Axis(fig_out[1,1], title=ax.title[], xlabel=ax.xlabel[], ylabel=ax.ylabel[])
    for obs in traj_observables
        obs(ax_out)   # caller provides a lambda that re-plots
    end
    mkpath(dirname(fname))
    save(fname, fig_out, px_per_unit=2)
    @info "Saved $fname"
end


# --- simulation variables ---
E  = Observable(Float64(init_var[1]))
x0 = Observable(Float64(init_var[2]))
y0 = Observable(Float64(init_var[3]))
py0= Observable(Float64(init_var[4]))
T  = Observable(  Int64( num_int[1]))

# --- analysis observables ---

y_sec_obs  = Observable(Float64[]) # surface of section
py_sec_obs = Observable(Float64[])

y0_range_obs  = Observable(Float64[])
py0_range_obs = Observable(Float64[])

# --- figure parameters---
ideal_size= (1920 , 1200)

with_theme(QUARTO_THEME) do
    init = lift(E, x0, y0, py0, T) do e, x, y, py, T
        a, m, w = param
        v       = potential(x, y; p=param)  # potential energy
        psquare = 2m*(e - v)  # p_y^2 if p_x=0 
        pmax    = sqrt(max(0.0, psquare))  # results in a maximal p_y
        py_c    = clamp(py, 0.0, pmax)  # clamp the observable py
        pdiff   = psquare - py_c^2  
        px      = sqrt(max(0.0, pdiff))
        # println(" – the polynomial\nsolve the polynomial for y0min, y0max:\n$(limit_of_initial_y0(e))")
        (x, y, px, py_c, T)
    end

    # --- derived trajectory ---
    sol = lift(init) do (x, y, px, py, T)
        y_sec, py_sec = Float64[], Float64[]
        cb  = section_callback(y_sec, py_sec)
        s   = solve_trajectory(
                [x, y, px, py], 
                tspan=(0.0, T), 
                dt=num_int[2], 
                p=param,
                callback=cb
                )
        y_sec_obs[] =   y_sec
        py_sec_obs[]=   py_sec
        s
    end

    # --- (phase space) trajectory ---
    traj_x = lift(s -> s[1, :], sol)   # all x values: row 1, all timesteps
    traj_y = lift(s -> s[2, :], sol)   # all y values: row 2, all timesteps


    # --- layout ---
    fig = Figure(size=ideal_size)
    # real space trajectory
    ax1 = Axis(fig[1:2,1], title="Space (x, y)", xlabel="x", ylabel="y")
    xlims!(ax1, -1.0, 1.0)
    ylims!(ax1, -1.0, 1.0)
    lines!(ax1, traj_x, traj_y, linewidth=1.0)
    # plotting the contour of potential into fig 1
    r = range(-1.0, 1.0, length = 120)
    levels = logrange(5.0*0.009,6.9*0.089,7)
    epot = [potential(x,y; p=param) for x in r, y in r]
    contour!(ax1, r, r, epot, labels=true, levels=levels, colormap=:hsv, colorscale=identity)

    # phase space trajectory
    ax2 = Axis(fig[1:2,2], title="Section (x=0)", xlabel="y", ylabel="p_y")
    scatter!(ax2, y_sec_obs, py_sec_obs, markersize=5)
    # calculating maximum cricumfarence of phase space trajectory

    on(E; update=true) do e
        y0_roots        = real.(filter(r -> abs(imag(r)) < 1e-10, limit_of_initial_y0(e)))
        y0_1, y0_2 = extrema(y0_roots)  # y0_roots[1], y0_roots[2]
        y0_range_obs[]  = collect(range(y0_1, y0_2, length=100))
        py0_range_obs[] = limit_of_initial_py0(y0_range_obs[], e)
        _, p_max    = extrema(py0_range_obs[])
        xlims!(ax2, y0_1,y0_2)
        ylims!(ax2, -p_max, p_max)
    end
    scatter!(ax2, y0_range_obs, py0_range_obs, markersize=5, color=C_CREAM)
    scatter!(ax2, y0_range_obs, lift(v -> -v, py0_range_obs), markersize=5, color=C_CREAM)


    # x(t)
    sol_t = lift(s->s.t, sol)
    ax3 = Axis(fig[3,1:2], title="x(t)", xlabel="t", ylabel="x")
    lines!(ax3, sol_t, traj_x, linewidth=0.5)
    ylims!(ax3,  -1.0, 1.0)
    on(T; update=true) do t
        xlims!(ax3, 0.0, t)
    end

    # y(t)
    ax4 = Axis(fig[4,1:2], title="y(t)", xlabel="t", ylabel="y")
    lines!(ax4, sol_t, traj_y, linewidth=0.5)
    ylims!(ax4,  -1.0, 1.0)
    on(T; update=true) do t
        xlims!(ax4, 0.0, t)
    end

    # --- controls ---
    sg = SliderGrid(fig[1, 3],
        (label="Energy E",  range=0.01:0.0001:0.3,   startvalue=0.083),
        (label="x₀",        range=-1.0:0.01:1.0,    startvalue=0.0),
        (label="y₀",        range=-0.5:0.01:1.0,    startvalue=0.0),
        (label="py₀",       range=0.0:0.01:1.0,     startvalue=0.0),
        (label="T",       range=100:50:2000,     startvalue=num_int[1])
    )
    on(sg.sliders[1].value) do v; E[]   = v; end
    on(sg.sliders[2].value) do v; x0[]  = v; end
    on(sg.sliders[3].value) do v; y0[]  = v; end
    on(sg.sliders[4].value) do v; py0[] = v; end
    on(sg.sliders[5].value) do v; T[]   = v; end

    # --- save figure/data ---
    # Label(fig[2, 3][1,1:4], "save figure")
    # save_fig_btn_traj = Button(fig[2, 3][1,1], label="traj")
    # save_fig_btn_phasesect = Button(fig[2, 3][1,2], label="section")
    # save_fig_btn_xt = Button(fig[2, 3][1,3], label="x(t)")
    # save_fig_btn_yt = Button(fig[2, 3][1,4], label="y(t)")

    # on(save_fig_btn.clicks) do _
    #     fname = joinpath(@__DIR__, "../../figures/henon-heiles", "explore_E$(round(E[], digits=4)).png")
    #     mkpath(dirname(fname))
    #     # save only ax1 (trajectory plot)
    #     fig_export = Figure(size=ideal_size)
    #     ax_export = Axis(fig_export[1,1])
    #     lines!(ax_export, traj_x[], traj_y[], linewidth=1.0)
    #     save(fname, fig_export, px_per_unit=2)
    # end




    Label(fig[5, 1:2], lift(v -> "effective py₀ = $(round(v[4], digits=4))", init))

    display(GLMakie.Screen(), fig)
end