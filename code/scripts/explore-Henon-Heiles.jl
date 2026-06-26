import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using GLMakie
using JLD2

include("../models/henon_heiles.jl")
using .HenonHeiles
include("../styles/makie_theme.jl")

# checking for paths
mkpath(joinpath(@__DIR__, "../../figures/henon-heiles/explore"))
mkpath(joinpath(@__DIR__, "../../data/henon-heiles/explore"))

# --- simulation variables ---
cfg   = load_config(joinpath(@__DIR__, "../sim_config/henon_heiles.json"))
param = (cfg.a, cfg.m, cfg.w)


E  = Observable(Float64(cfg.E))
x0 = Observable(Float64(cfg.x0))
y0 = Observable(Float64(cfg.y0))
py0= Observable(Float64(cfg.py0))
T  = Observable(  Int64(cfg.T))

# --- analysis observables ---
y_sec_obs  = Observable(Float64[]) # surface of section
py_sec_obs = Observable(Float64[])

y0_range_obs  = Observable(Float64[])
py0_range_obs = Observable(Float64[])

# --- figure parameters---
full_screen_size= (1920 , 1200)
save_size       = (1200, 1200)
# --- contour plot ---
r      = range(-1.0, 1.0, length=120)
levels = logrange(5.0*0.009, 6.9*0.089, 7)
epot   = [potential(x,y, param) for x in r, y in r]


# --- functions ---

function bind_textbox(fig_pos, slider, label_text; parser=Float64)
    tb = Textbox(fig_pos, placeholder=label_text,
                 validator = s -> tryparse(parser, s) !== nothing, width=20)
    on(tb.stored_string) do s
        v = tryparse(parser, s)
        v === nothing && return
        set_close_to!(slider, v)
    end
    tb
end

# === User interface ===
with_theme(QUARTO_THEME) do
    init = lift(E, x0, y0, py0, T) do e, x, y, py, T
        a, m, w = param
        v       = potential(x, y, param)  # potential energy
        psquare = 2m*(e - v)  # p_y^2 if p_x=0 
        pmax    = sqrt(max(0.0, psquare))  # results in a maximal p_y
        py_c    = clamp(py, 0.0, pmax)  # clamp the observable py
        pdiff   = psquare - py_c^2  
        px      = sqrt(max(0.0, pdiff))
        (x, y, px, py_c, T)
    end

    # --- derived trajectory ---
    sol = lift(init) do (x, y, px, py, T)
        y_sec, py_sec = Float64[], Float64[]
        cb  = section_callback(y_sec, py_sec)
        s   = solve_trajectory(
                [x, y, px, py], 
                param,
                (0.0, T), 
                cfg.dt; 
                callback=cb
                )
        y_sec_obs[] =   y_sec
        py_sec_obs[]=   py_sec
        s
    end

    param_label = lift(E, init) do e, (x, y, px, py, T)
        a, m, w = param
        """
        a  = $a
        m  = $m
        ω  = $w
        T  = $T
        dt = $(cfg.dt)
        E  = $e
        x₀ = $x
        y₀ = $y
        py₀= $py
        px₀= $px
        """
    end



    # --- (phase space) trajectory ---
    traj_x = lift(s -> s[1, :], sol)   # all x values: row 1, all timesteps
    traj_y = lift(s -> s[2, :], sol)   # all y values: row 2, all timesteps


    # === layout ===
    fig = Figure(size=full_screen_size)

    # == Plot: configuration space ==
    ax1 = Axis(fig[1:2,1], title="Configuration space", xlabel="x", ylabel="y")
    xlims!(ax1, -1.0, 1.0)
    ylims!(ax1, -0.6, 1.2)
    lines!(ax1, traj_x, traj_y, linewidth=1.0)

    contour!(ax1, r, r, epot, labels=true, levels=levels, colormap=:hsv, colorscale=identity)

    Label(fig[2, 3][1,1:4], "Save options")
    # --- save trajectory figure ---
    save_fig_btn_traj = Button(fig[2, 3][2,1], label="traj", labelcolor=:black)
    on(save_fig_btn_traj.clicks) do _
        println("--- save trajecctory plot ---")
        println("parameters:\n$(param_label[])")
        fname = joinpath(@__DIR__, "../../figures/henon-heiles/explore", "traj-$(round(E[], digits=4))-$(x0[])-$(y0[])-$(py0[]).png")
        with_theme(QUARTO_THEME) do
            fig_export = Figure(size=save_size)
            ax_export  = Axis(fig_export[1,1])
            lines!(ax_export, traj_x[], traj_y[], linewidth=1.0)
            contour!(ax_export, r, r, epot, labels=true, levels=levels, colormap=:hsv, colorscale=identity)
            xlims!(ax_export, -1, 1)
            ylims!(ax_export, -0.6, 1)
            ax_export.xlabel = "x"
            ax_export.ylabel = "y"
            save(fname, fig_export, px_per_unit=2)
        end
        println("--- done: $fname ---")
    end

    # == Plot: phase space trajectory ==
    ax2 = Axis(fig[1:2,2], title="Section (x=0)", xlabel="y", ylabel="p_y")
    scatter!(ax2, y_sec_obs, py_sec_obs, markersize=5)  # from trajectory
    # calculating maximum cricumfarence of phase space trajectory
    on(E; update=true) do e
        y0_roots        = real.(filter(r -> abs(imag(r)) < 1e-10, limit_of_initial_y0(e, param)))
        y0_1, y0_2 = y0_roots[1], y0_roots[2]
        dy = (y0_roots[2] - y0_roots[1])/100
        y0_range       = collect(range(y0_1, y0_2, length=100))
        y0_range_obs[] = y0_range
        py0_range_obs[] = limit_of_initial_py0(y0_range, e, param)
        _, p_max    = extrema(py0_range_obs[])
        dp = p_max/50
        xlims!(ax2, y0_1-dy, y0_2+dy)
        ylims!(ax2, -p_max-dp, p_max+dp)
    end
    scatter!(ax2, y0_range_obs, py0_range_obs, markersize=5, color=C_CREAM)  # max py
    scatter!(ax2, y0_range_obs, lift(v -> -v, py0_range_obs), markersize=5, color=C_CREAM) # min -py


    # --- save surface section figure ---
    save_fig_btn_phasesect = Button(fig[2, 3][2,2], label="section", labelcolor=:black)
    on(save_fig_btn_phasesect.clicks) do _
        println("--- save surface of section plot ---")
        println("parameters:\n$(param_label[])")
        fname = joinpath(@__DIR__, "../../figures/henon-heiles/explore", "surf_of_sec-$(round(E[], digits=4))-$(x0[])-$(y0[])-$(py0[]).png")
        with_theme(QUARTO_THEME) do
            fig_export = Figure(size=save_size)
            ax_export  = Axis(fig_export[1,1])
            scatter!(ax_export, y_sec_obs[], py_sec_obs[], markersize=5)
            scatter!(ax_export, y0_range_obs[], py0_range_obs[], markersize=5, color=C_CREAM)
            scatter!(ax_export, y0_range_obs[], -py0_range_obs[], markersize=5, color=C_CREAM)
            y_min, y_max = extrema(y0_range_obs[])
            p_max        = maximum(py0_range_obs[])
            dy = (y_max - y_min) / 50
            dp = p_max / 30
            xlims!(ax_export, y_min-dy, y_max+dy)
            ylims!(ax_export, -p_max-dp, p_max+dp)
            ax_export.xlabel = "y"
            ax_export.ylabel = "p_y"
            save(fname, fig_export, px_per_unit=2)
        end
        println("--- done: $fname ---")
    end

    # == Plot: x(t) ==
    sol_t = lift(s->s.t, sol)
    ax3 = Axis(fig[3,1:2], title="x(t)", xlabel="t", ylabel="x")
    lines!(ax3, sol_t, traj_x, linewidth=0.5)
    ylims!(ax3,  -1.0, 1.0)


    # --- save x(t) figure ---
    save_fig_btn_xt = Button(fig[2, 3][2,3], label="x(t)", labelcolor=:black)
    on(save_fig_btn_xt.clicks) do _
        println("--- save xt plot ---")
        println("parameters:\n$(param_label[])")
        fname = joinpath(@__DIR__, "../../figures/henon-heiles/explore", "xt-$(round(E[], digits=4))-$(x0[])-$(y0[])-$(py0[]).png")
        with_theme(QUARTO_THEME) do
            fig_export = Figure(size=save_size)
            ax_export  = Axis(fig_export[1,1])
            lines!(ax_export, sol_t[], traj_x[], linewidth=0.5)
            ylims!(ax_export, -1.0, 1.0)
            xlims!(ax_export, 0.0, T[]+0.01)
            ax_export.xlabel = "t"
            ax_export.ylabel = "x(t)"
            save(fname, fig_export, px_per_unit=2)
        end
        println("--- done: $fname ---")
    end

    # == Plot: y(t) ==
    ax4 = Axis(fig[4,1:2], title="y(t)", xlabel="t", ylabel="y")
    lines!(ax4, sol_t, traj_y, linewidth=0.5)
    ylims!(ax4,  -1.0, 1.0)
    # update the time axis of trajectories
    on(T; update=true) do t
        xlims!(ax3, 0.0, t+0.01)
        xlims!(ax4, 0.0, t+0.01)
    end


    # --- save y(t) figure ---
    save_fig_btn_yt = Button(fig[2, 3][2,4], label="y(t)", labelcolor=:black)
    on(save_fig_btn_yt.clicks) do _
        println("--- save yt plot ---")
        println("parameters:\n$(param_label[])")
        fname = joinpath(@__DIR__, "../../figures/henon-heiles/explore", "yt-$(round(E[], digits=4))-$(x0[])-$(y0[])-$(py0[]).png")
        with_theme(QUARTO_THEME) do
            fig_export = Figure(size=save_size)
            ax_export  = Axis(fig_export[1,1])
            lines!(ax_export, sol_t[], traj_y[], linewidth=0.5)
            ylims!(ax_export, -1.0, 1.0)
            xlims!(ax_export, 0.0, T[]+0.01)
            ax_export.xlabel = "t"
            ax_export.ylabel = "y(t)"
            save(fname, fig_export, px_per_unit=2)
        end
        println("--- done: $fname ---")
    end

    # === controls of initial parameters ===
    sg = SliderGrid(fig[1, 3],
        (label="E",  range=0.01:0.0001:0.16667,   startvalue=0.083),
        (label="x₀",        range=-1.0:0.01:1.0,    startvalue=0.0),
        (label="y₀",        range=-0.5:0.01:1.0,    startvalue=0.0),
        (label="py₀",       range=0.0:0.01:1.0,     startvalue=0.0),
        (label="T",       range=500:500:100000,     startvalue=5000)
    )

    on(sg.sliders[1].value) do v; E[]   = v; end
    on(sg.sliders[2].value) do v; x0[]  = v; end
    on(sg.sliders[3].value) do v; y0[]  = v; end
    on(sg.sliders[4].value) do v; py0[] = v; end
    on(sg.sliders[5].value) do v; T[]   = v; end


    # tb_E   = bind_textbox(sg.layout[1, 4], sg.sliders[1], "E")
    # tb_x0  = bind_textbox(sg.layout[2, 4], sg.sliders[2], "x₀")
    # tb_y0  = bind_textbox(sg.layout[3, 4], sg.sliders[3], "y₀")
    # tb_py0 = bind_textbox(sg.layout[4, 4], sg.sliders[4], "py₀")
    # tb_T   = bind_textbox(sg.layout[5, 4], sg.sliders[5], "T"; parser=Int64)
    
    # on(E)   do v; tb_E.displayed_string[]   = string(round(v, digits=4)); end
    # on(x0)  do v; tb_x0.displayed_string[]  = string(v); end
    # on(y0)  do v; tb_y0.displayed_string[]  = string(v); end
    # on(py0) do v; tb_py0.displayed_string[] = string(v); end
    # on(T)   do v; tb_T.displayed_string[]   = string(v); end

    # displaying all parameters
    Label(fig[3, 3], param_label, justification=:left, font="JetBrainsMono Nerd Font")
    is_valid = lift(init) do (x, y, px, py, T)
        px == 0.0 && py == 0.0 ? "⚠ outside energy surface" : ""
    end
    Label(fig[5, 1:2], is_valid, color=:red)

    # --- save data ---
    save_data_bin_btn = Button(fig[2, 3][3,1:4], label="save data", labelcolor=:black)
    on(save_data_bin_btn.clicks) do _
        println("--- save data ---")
        println("parameters:\n$(param_label[])")
        fname1 = joinpath(@__DIR__, "../../data/henon-heiles/explore",
            "sim_traj-T$(T[])-E$(round(E[], digits=4))-x0$(x0[])-y0$(y0[])-py0$(py0[]).jld2")
        fname2 = joinpath(@__DIR__, "../../data/henon-heiles/explore",   
            "surf_of_sec-T$(T[])-E$(round(E[], digits=4))-x0$(x0[])-y0$(y0[])-py0$(py0[]).jld2")
        jldsave(fname1; t=sol_t[], u=sol[].u)  # saving trajectories
        jldsave(fname2; y=y_sec_obs[], py=py_sec_obs[])  # saving surface of section
        println("--- done: $fname1 and $fname2 ---")
    end



    display(GLMakie.Screen(), fig)
end