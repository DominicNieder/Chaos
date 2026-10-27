const FIG_DIR    = joinpath(@__DIR__, "../../figures/two-boxed-charges/")



"""
takes a bunch of initial conditions u0s and returns three plots: the trajectories in configuration space, the surface of section, and the energy vs time.

    u0s: a vector of initial conditions, each of which is a 4-vector [x1, x2, p1, p2]
    p: the parameters of the system, a NamedTuple with fields C, m1, m2, L1, L2, del
    E: the energy of the system
    int_time: the integration time for the trajectories
    style: the style of the plots, either :print or :presentation
    save_fig: whether to save the figures to disk   
"""
function explore_1D_boxes(seeds, p, E; int_time=10_000, style=:print, save_fig=false)
    u0s = [lift(v,E,p) for v in seeds]
    set_style!(style)
    
    f = Figure(size=(1300,900))
    ax_sec    = Axis(f[2,1:2], xlabel=L"x_2\, [L]", ylabel=L"p_2")
    ax_traj   = Axis(f[1,1:2], xlabel=L"\hat{x}_1\, [\text{L}]", ylabel=L"\hat{x}_2 \, [\text{L}]")
    ax_energy = Axis(f[2,3],   xlabel=L"t", ylabel=L"\text{Energy}\,[\mathcal{C}/L]")
    # [hlines!(ax_traj, get_boxes(p)[j][i], color=INK, linestyle=:dash) for i in 1:2, j in 1:2]

    Es0 = Float64[]                                  # initial energy per seed, for the legend
    for (i, u0) in enumerate(u0s)
        E0 = energy(u0, p); push!(Es0, E0)
        u, t, pts = get_traj(u0, int_time; p=p, abstol=INT_TOL, reltol=INT_TOL)

        xs   = first.(u, 2)
        xone = Point2f.(zip(t, first.(xs)))
        xtwo = Point2f.(zip(t, last.(xs)))
        sec  = Point2f.(first.(pts, 2))

        lines!(ax_traj, Point2f.(xs), color=(pick_color(i),0.8))
        # lines!(ax_traj, xtwo[1:2000], color=pick_color(i), linestyle=:dash)
        scatter!(ax_sec, sec, color=pick_color(i), markersize=4)
        lines!(ax_energy, t, [abs(energy(ui,p) - E)*p.L1/E/p.C for ui in u], color=pick_color(i))
    end
    scatter!(ax_sec, boundary(E; p,n=1000), color=INK, markersize=3)

    # ── legend: component (line style) + seeds (colour, with E₀) ──
    comp_elems = [LineElement(color=:gray70, linestyle=:solid, linewidth=2),
                LineElement(color=:gray70, linestyle=:dash,  linewidth=2),
                LineElement(color=INK,     linestyle=:dash,  linewidth=1.5),
                MarkerElement(color=INK, marker=:circle, markersize=6)]
    comp_labels = [L"x_1", L"x_2", "box edges", "section boundary"]

    seed_elems  = [LineElement(color=pick_color(i), linewidth=3) for i in eachindex(u0s)]
    seed_labels = ["seed $(seeds[i])   " for i in eachindex(u0s)]

    Legend(f[1,3], [comp_elems, seed_elems], [comp_labels, seed_labels],
        ["Component", "Seeds  (E = $(round(E, sigdigits=4)))"];
        framevisible   = false,
        tellheight     = false,
        halign         = :left, valign = :top,
        titlehalign    = :left, titlesize = 12,
        gridshalign    = :left,
        labelsize      = 11,
        patchsize      = (18, 10),
        rowgap         = 1,
        groupgap       = 8,
        padding        = (2, 2, 2, 2))

    colsize!(f.layout, 3, Auto(0.45))   # keep column 3 narrow so the legend sits snug
    ylims!(ax_sec, -200,200)
    folder = joinpath(FIG_DIR, "explore_boxGap/")
    if save_fig
        save(joinpath(folder, "E$E-P$p.png"), f)   # save AFTER the legend exists
    end
    return (;f, ax_sec, ax_traj, ax_energy)
end


function explore_del_between_boxes(ps; int_time=10_000, style=:print, save_fig=false)
    seeds = [[1.54e-4, 20], [0.0006,0.0], [5.3e-4,18.0], [0.00015,40],  [1e-4,0.0],  [7.7735e-4, 0.0],  [6e-5, 63.0]]  # C<0
    out = []
    for (i, p_i) in enumerate(ps)
        E  = 500 / (p_i.C/p_i.L1)
        
        fi = explore_1D_boxes(seeds, p_i, E; int_time=int_time, style=style, save_fig=save_fig)
        push!(out, fi)
    end
    return out
end



function explore_energies(seeds, E0; p=(;C=-1, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= 1e-15), int_time=10_000, style=:print, save_fig=false)
    out = []
    # here I want to follow the seeds that I find at E0
    for E in 10 .^ range(-3, stop=0, length=10)
        fi = explore_1D_boxes(seeds, p, E; int_time=int_time, style=style, save_fig=save_fig)
        push!(out, fi)
    end
    return out
end

"""
    indent_jump(theta; del=0.1)

Return a Boolean vector matching `theta`. For each adjacent pair, mark the
point after it when `abs(theta[i] - theta[i + 1]) > del`. Non-finite
differences are not marked.
"""
function indent_jump(θ; del=0.05)
    isfinite(del) && del >= 0 || throw(ArgumentError("del must be finite and nonnegative"))
    jumps = falses(length(θ))
    for i in 1:(length(θ) - 1)
        difference = abs(θ[i] - θ[i + 1])
        isfinite(difference) && difference > del && (jumps[i + 1] = true)
    end
    return jumps
end

"""
    analyse_followed_traj(orbs, p; size=(1300, 900), trajectory_every=100, style=:dark)

Compute monodromy diagnostics for a set of followed orbits and visualise their
continuation in phase and configuration space.

The function evaluates the monodromy matrix for each orbit, computes the
associated phase variable `theta`, detects jumps in this phase, and prints a
summary of the largest discontinuity. It then generates a continuation plot and
a configuration-space plot, highlighting orbit points with phase jumps.

Returns a named tuple:
- `f`, `ax`: phase-space continuation figure and axis
- `config_f`, `config_ax`: configuration-space figure and axis
- `orbits`: original orbit data augmented with `M`, `check`, `theta`, and `jump`
- `monodrome`: monodromy/eigenvalue plots
"""
function analyse_followed_traj(orbs, p; 
       size=(1300, 900), 
       trajectory_every=100, 
       style=:dark)

    trajectory_every > 0 || throw(ArgumentError("trajectory_every must be positive"))
    ana = monodrome(orbs; p=p, verbose=false)
    orbits_with_mono = copy(orbs)
    orbits_with_mono.M = ana.M
    orbits_with_mono.check = ana.check
    orbits_with_mono.theta = [
        valid ? get_phase2(eigvals(M)) / (2π) : NaN
        for (M, valid) in zip(ana.M, ana.check)
    ]
    orbits_with_mono.jump = indent_jump(orbits_with_mono.theta; del=0.005)
    njumps = count(orbits_with_mono.jump)
    theta_differences = abs.(diff(orbits_with_mono.theta))
    finite_difference_indices = findall(isfinite, theta_differences)
    println(any(orbits_with_mono.jump) ?
            "Detected $njumps jump(s)." :
            "No jumps detected.")
    if isempty(finite_difference_indices)
        println("No finite adjacent theta differences to inspect.")
    else
        max_difference_index = finite_difference_indices[
            argmax(theta_differences[finite_difference_indices])
        ]
        println("Largest adjacent theta change: ",
                theta_differences[max_difference_index],
                " between E=", orbits_with_mono.E[max_difference_index],
                " and E=", orbits_with_mono.E[max_difference_index + 1],
                " (jump threshold del=0.05)")
    end
    
    eigenvalue_plots = graphs(orbits_with_mono;
              fsize=size, 
              style=style, 
              save_fig=false)

    order = sortperm(orbs.E)
    points = [Point2f(orbs.v[i][1], orbs.v[i][2]) for i in order]
    jump_points = orbits_with_mono.jump[order]
    length(points) >= 2 || error("need at least two orbit points to plot the energy flow")

    f = Figure(size=size)
    ax = Axis(f[1, 1], xlabel=L"x_2", ylabel=L"p_2",
              title="Orbit continuation (arrows show increasing energy)")
    scatter!(ax, points[.!jump_points]; color=pick_color(1), markersize=3)

    arrow_indices = 1:max(1, cld(length(points) - 1, 30)):(length(points) - 1)
    arrow_positions = points[arrow_indices]
    arrow_directions = [points[i + 1] - points[i] for i in arrow_indices]
    arrows2d!(ax, arrow_positions, arrow_directions;
              color=pick_color(1), tipwidth=12, tiplength=12, shaftwidth=1.5)
    scatter!(ax, points[jump_points]; color=INK, markersize=9,
             strokecolor=pick_color(1), strokewidth=1)

    config_f = Figure(size=size)
    config_ax = Axis(config_f[1, 1], xlabel=L"x_1", ylabel=L"x_2",
                     aspect=DataAspect(),
                     title="Configuration-space orbits (first and every $trajectory_every)")
    regular_indices = 1:trajectory_every:length(order)
    jump_indices = order[findall(jump_points)]
    selected_orbit_indices = sort!(unique(vcat(order[regular_indices], jump_indices)))
    for i in selected_orbit_indices
        u0 = lift(orbs.v[i], orbs.E[i], p)
        u0 === nothing && error("could not lift orbit point at E=$(orbs.E[i])")
        trajectory, _, _ = get_traj(u0, orbs.T[i]; p)
        is_jump = orbits_with_mono.jump[i]
        lines!(config_ax, Point2f.(getindex.(trajectory, 1), getindex.(trajectory, 2));
               color=is_jump ? INK : (pick_color(1), 0.6),
               linewidth=is_jump ? 3 : 2)
    end
    if !isempty(jump_indices)
        jump_states = [lift(orbs.v[i], orbs.E[i], p) for i in jump_indices]
        any(isnothing, jump_states) && error("could not lift a jump orbit point")
        scatter!(config_ax, [Point2f(u[1], u[2]) for u in jump_states];
                 color=INK, markersize=10, strokecolor=pick_color(1), strokewidth=1)
    end
    return (;f, ax, config_f, config_ax, orbits=orbits_with_mono,
            monodrome=eigenvalue_plots)
end
