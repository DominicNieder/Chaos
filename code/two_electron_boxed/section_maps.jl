using Roots

uc(c) = find_zero(u -> 2/π*(atan(u) + u/(1+u^2)) - c, (0.0, 1e6), Roots.Bisection())
ε(E, m, c) = uc(c) * sqrt(-2m*E)


"""
    section_map(E, p; nx=5, np=5, tint=1_000, plims=nothing)

Seeds the Poincaré section at energy `E` (model units), integrates each seed
for `tint`, and returns:
  - `seed.f`  : seeds on the section with the energy boundary
  - `sec.f`   : surface of section (left), energy drift (right top), legend (right bottom)
  - raw data  : `seeds` (kept), `sec_map`
"""
function section_map(E, p; nx=5, np=5, tint=1_000, plims=(-ε(E, 1, 0.99),ε(E, 1, 0.99)))
    seeds     = sample_grid(E, p; nx, np)
    sec_bound = boundary(E; p, n=500)

    f    = Figure(size=(1300, 900))
    ax   = Axis(f[1:2, 1:2], xlabel=L"x_2\,[L]", ylabel=L"p_2",
                title="Surface of section, E = $(round(E; sigdigits=4))")
    axen = Axis(f[1, 3], xlabel=L"t", ylabel=L"|E(t)-E|\,/\,|E|",
                yscale=log10, title="Energy drift")
    Escale = E == 0 ? abs(p.C / p.L1) : abs(E)
    scatter!(ax, sec_bound, color=INK, markersize=3)

    # ── integrate ──
    kept    = Vector{Float64}[]
    sec_map = Vector{Point2f}[]

    naborted = 0
    for seed in seeds
        u0 = lift(seed, E, p)
        u0 === nothing && continue
        u, ts, pts = get_traj(u0, tint; p)
        # if ts[end] < tint                      # solver aborted (near-collision, dt → eps)
        #     naborted += 1
        #     continue
        # end
        push!(kept, seed)
        sec = Point2f.(first.(pts, 2))
        push!(sec_map, sec)
        scatter!(ax, sec, color=pick_color(length(kept)), markersize=4)
        err = [max(abs(energy(ui, p) - E) / Escale, eps()) for ui in u]
        lines!(axen, ts, err, color=pick_color(length(kept)))
    end

    naborted > 0 && @info "E = $E: $naborted of $(length(seeds)) seeds aborted by the solver"
    isempty(kept) && error("all seeds aborted at E = $E")

    # ── seed figure ──
    fseed    = Figure(size=(1300, 900))
    ax_seeds = Axis(fseed[1, 1], xlabel=L"x_2\,[L]", ylabel=L"p_2",
                    title="Seeds, E = $(round(E; sigdigits=4))")
    scatter!(ax_seeds, sec_bound, color=INK, markersize=3)
    scatter!(ax_seeds, Point2f.(kept), color=[pick_color(i) for i in eachindex(kept)],
             markersize=10)
    plims === nothing || ylims!(ax_seeds, plims...)

    # ── section (left) + energy drift and legend (right) ──
    plims === nothing || ylims!(ax, plims...)

    # ── legend ──
    comp_elems  = [MarkerElement(color=INK, marker=:circle, markersize=6)]
    seed_elems  = [MarkerElement(color=pick_color(i), marker=:circle, markersize=8)
                   for i in eachindex(kept)]
    seed_labels = [@sprintf("(%.3g, %.3g)", s[1], s[2]) for s in kept]
    Legend(f[2, 3], [comp_elems, seed_elems], [["section boundary"], seed_labels],
           ["Component", "Seeds (x₂, p₂)"];
           framevisible=false, tellheight=false,
           halign=:left, valign=:top, titlehalign=:left, titlesize=12,
           gridshalign=:left, labelsize=11, rowgap=1, groupgap=8,
           nbanks = length(kept) > 12 ? 2 : 1)
    colsize!(f.layout, 3, Auto(0.45))

    return (; E, seed=(; f=fseed, ax=ax_seeds), sec=(; f, ax, axen),
              seeds=kept, sec_map)
end


"""
For attractive section maps; I am expecting the system not to change, I expect the system to be topologicaly the same. I expect the monodromy method to result be constand M(E)= M
"""


"""
    scan_section_maps(Es, p; nx=5, np=5, tint=1_000, plims=nothing,
                      save_fig=false, folder=joinpath(FIG_DIR, "section_maps"))

Runs `section_map` for every energy in `Es` (model units). Energies without an
accessible section are skipped with a warning. Returns a vector of results;
view with `display(res[i].sec.f)` or `display(res[i].seed.f)`.
"""
function scan_section_maps(Es, p; nx=5, np=5, tint=1_000,
                           save_fig=false,
                           folder=joinpath(FIG_DIR, "section_maps"),c=0.99)
    save_fig && mkpath(folder)
    tag = @sprintf("C%+g_del%.0e_L%g-%g", p.C, p.del, p.L1, p.L2)

    out = []
    @showprogress desc="section maps" for E in Es
        res = try
            plims = (-ε(E, 1.0, c), ε(E, 1.0, c))
            section_map(E, p; nx, np, tint, plims)
        catch err
            @warn "skipped E = $E" exception=(err, catch_backtrace())
            continue
        end
        push!(out, res)

        if save_fig
            Etag = @sprintf("E%+.4g", E)
            save(joinpath(folder, "seeds_$(Etag)_$(tag).png"),   res.seed.f)
            save(joinpath(folder, "section_$(Etag)_$(tag).png"), res.sec.f)
        end
    end
    return out
end



# ============================================================================================
#        Periodic-orbit search up to period nmax — reuses the Hénon–Heiles pipeline:
#        sample_grid → sweep!/analyse_seed → minPeriodicity → already_found → monodrome
# ============================================================================================

"`orbit_table()` plus the columns `already_found` expects (prime, sec_y, sec_py)."
function po_table()
    df = orbit_table()
    df.prime  = Int[]
    df.sec_y  = Vector{Float64}[]
    df.sec_py = Vector{Float64}[]
    return df
end

kind_marker(k) = k == 0 ? :star5 : KIND_MS[k]
kind_label(k)  = k == 0 ? "monodromy failed" : KIND_LABEL[k]


"""
    find_periodic_orbits(E, p; nmax=8, nx=5, np=5, ...)

For n = 1…nmax: roots of Tⁿ(v) = v via `sweep!` (i.e. `analyse_seed`),
prime period and section points via `minPeriodicity`, duplicates removed with
`already_found`. Stability from `monodrome` (4×4); transverse multipliers are extracted with
`get_eigenvals` and classified with `monodromy_kind_index`.

Columns: E, v, T, str ("N<period>"), prime, sec_y, sec_py, M, check, τ, kind, ρ, lyap
"""
function find_periodic_orbits(E, p; nmax=8, nx=5, np=5, tmax=100_000.0,
                              dup_tol=1e-7, verbose=false, show_progress=true)
    seeds = sample_grid(E, p; nx, np)
    prm   = SectionParams(E, p; tmax, nfast=1, ndense=40,
                          save_everystep=false, save_start=false)
    df    = po_table()

    for n in 1:nmax
        prm.nmax_fast[] = n                              # stop after n crossings
        roots = sweep!(orbit_table(), seeds, n, prm; verbose, show_progress)
        for o in eachrow(roots)
            mp = minPeriodicity(o.v, prm; search=max(40, 2nmax))
            mp.Nperiod === nothing && continue
            already_found(df, o.v, mp.Nperiod; pmap_prime_tol=dup_tol) && continue
            push!(df, (; E=o.E, v=o.v, T=o.T, str="N$(mp.Nperiod)", prime=mp.Nperiod,
                         sec_y=mp.pMap[1, :], sec_py=mp.pMap[2, :]))
        end
    end
    isempty(df) && return df

    # ── stability (existing monodromy code) ──
    mono      = monodrome(df; p, verbose, show_progress)
    df.M      = mono.M
    df.check  = mono.check
    
    df.τ      = [c ? tr(M) - 2 : NaN for (M, c) in zip(df.M, df.check)]
    df.kind   = [c ? monodromy_kind_index(M) : 0 for (M, c) in zip(df.M, df.check)]
    df.ρ      = [k == 1 ? abs(angle(get_eigenvals(M)[1])) / 2π : NaN
                 for (M, k) in zip(df.M, df.kind)]       # rotation number, elliptic only
    df.lyap   = [c ? get_lyapunov(M, T) : NaN for (M, T, c) in zip(df.M, df.T, df.check)]

    return sort!(df, [:prime, :T])
end

    


"""
    periodic_orbit_map(E, p; nmax=8, nx=5, np=5, plims=nothing,
                       background=true, nbg=3, tint=1_000, kwargs...)

Section with all periodic orbits (colour = prime period, marker = stability),
optionally on top of a `section_map` background, plus τ = tr(M) − 2 vs. T
(stable band |τ| < 2 shaded).
"""
function periodic_orbit_map(E, p; nmax=8, nx=5, np=5, plims=(-ε(E, 1, 0.99),ε(E, 1, 0.99)),
                            background=true, nbg=3, tint=1_000, orbit_df=nothing, kwargs...)
    df = orbit_df === nothing ? find_periodic_orbits(E, p; nmax, nx, np, kwargs...) : orbit_df
    background_points = Point2f[]

    f    = Figure(size=(1300, 900))
    ax   = Axis(f[1:2, 1:2], xlabel=L"x_2\,[L]", ylabel=L"p_2",
                title="Periodic orbits (N ≤ $nmax), E = $(round(E; sigdigits=4))")
    axtr = Axis(f[2, 3], xlabel=L"T", ylabel=L"\tau = \mathrm{tr}\,M - 2",
                yscale=Makie.pseudolog10, title="Stability")

    if background
        bg = section_map(E, p; nx=nbg, np=nbg, tint)
        for sec in bg.sec_map
            scatter!(ax, sec, color=(:gray60, 0.3), markersize=2)
            append!(background_points, sec)
        end
    end
    scatter!(ax, boundary(E; p, n=500), color=INK, markersize=3)
    hspan!(axtr, -2, 2, color=(:gray70, 0.3))

    for o in eachrow(df)
        c, m = pick_color(o.prime), kind_marker(o.kind)
        scatter!(ax, o.sec_y, o.sec_py, color=c, marker=m, markersize=10)
        o.check && scatter!(axtr, [o.T], [o.τ], color=c, marker=m, markersize=10)
    end
    plims === nothing || ylims!(ax, plims...)

    if isempty(df)
        text!(ax, 0.5, 0.5; text="no periodic orbits found", space=:relative,
              align=(:center, :center))
    else
        Ns = sort(unique(df.prime)); ks = sort(unique(df.kind))
        Legend(f[1, 3],
               [[MarkerElement(color=pick_color(N), marker=:circle, markersize=10) for N in Ns],
                [MarkerElement(color=:gray60, marker=kind_marker(k), markersize=10) for k in ks]],
               [["N = $N" for N in Ns], kind_label.(ks)],
               ["Period", "Stability"];
               framevisible=false, tellheight=false, halign=:left, valign=:top,
               titlehalign=:left, titlesize=12, gridshalign=:left, labelsize=11,
               rowgap=1, groupgap=8)
    end
    colsize!(f.layout, 3, Auto(0.45))

    return (; E, df, f, ax, axtr, background_points)
end


"""
    scan_periodic_orbits(Es, p; nmax=8, nx=5, np=5, plims=nothing,
                         save_fig=false, save_data=false, parallel=true,
                         show_progress=true, kwargs...)

Finds periodic orbits for each energy in parallel when multiple Julia threads
are available, then builds the Makie maps serially. Set `parallel=false` to
disable threaded searching. Progress is reported once per energy; set
`show_progress=false` to suppress it. Start Julia with `--threads=auto` to use
all available CPU cores.
Returns `(; maps, orbits)`: `display(maps[i].f)`; `orbits` holds all energies.
"""
function scan_periodic_orbits(Es, p; nmax=8, nx=5, np=5, plims=nothing,
                              save_fig=false, save_data=false,
                              background=true, nbg=3, tint=1_000,
                              folder=joinpath(FIG_DIR, "periodic_orbits"),
                              data_dir=DATA_DIR, parallel=true, show_progress=true, kwargs...)
    save_fig  && mkpath(folder)
    save_data && mkpath(data_dir)
    tag = @sprintf("C%+g_del%.0e_L%g-%g_N%d", p.C, p.del, p.L1, p.L2, nmax)

    energies = collect(Es)
    orbit_dfs = Vector{Union{Nothing,DataFrame}}(undef, length(energies))
    fill!(orbit_dfs, nothing)
    failures = Vector{Any}(undef, length(energies))
    fill!(failures, nothing)

    function scan_energy(i)
        E = energies[i]
        try
            orbit_dfs[i] = find_periodic_orbits(E, p; nmax, nx, np,
                                                show_progress=false, kwargs...)
        catch err
            failures[i] = (err, catch_backtrace())
        end
    end

    if parallel && Threads.nthreads() > 1 && length(energies) > 1
        @showprogress enabled=show_progress desc="energy scan" Threads.@threads for i in eachindex(energies)
            scan_energy(i)
        end
    else
        if parallel && length(energies) > 1 && Threads.nthreads() == 1
            @warn "Julia is running with one thread; start Julia with --threads=auto to parallelize the energy scan"
        end
        @showprogress enabled=show_progress desc="energy scan" for i in eachindex(energies)
            scan_energy(i)
        end
    end

    maps = []
    for (i, E) in enumerate(energies)
        if failures[i] !== nothing
            err, backtrace = failures[i]
            @warn "skipped E = $E" exception=(err, backtrace)
            continue
        end
        energy_plims = plims === nothing ? (-ε(E, 1, 0.99), ε(E, 1, 0.99)) : plims
        res = try
            periodic_orbit_map(E, p; nmax, nx, np, plims=energy_plims,
                               background, nbg, tint, orbit_df=orbit_dfs[i], kwargs...)
        catch err
            @warn "skipped plotting E = $E" exception=(err, catch_backtrace())
            continue
        end
        push!(maps, res)
        save_fig && save(joinpath(folder, @sprintf("po_E%+.4g_%s.png", E, tag)), res.f)
    end

    dfs    = [m.df for m in maps if !isempty(m.df)]
    orbits = isempty(dfs) ? po_table() : reduce(vcat, dfs; cols=:union)
    save_data && jldsave(joinpath(data_dir, "periodic_orbits_$(tag).jld2"); orbits, Es, p)
    return (; maps, orbits)
end


"""
    scan_periodic_orbits_interactive(Es, p; nmax=8, nx=5, np=5, kwargs...)

Run `scan_periodic_orbits` and display a GLMakie energy-slider view of the
periodic orbits on the section and their stability. Returns the scan maps,
combined orbit table `df` (with the `po_table()` columns), and the interactive
figure, axes, and slider.
"""
function scan_periodic_orbits_interactive(Es, p; nmax=8, nx=5, np=5,
                                          plims=nothing, background=false,
                                          save_fig=false, save_data=false,
                                          folder=joinpath(FIG_DIR, "periodic_orbits"),
                                          data_dir=DATA_DIR, kwargs...)
    set_style!(:dark)
    scan = scan_periodic_orbits(Es, p; nmax, nx, np, plims, save_fig, save_data,
                                folder, data_dir, background, kwargs...)
    isempty(scan.maps) && error("No energies produced a periodic-orbit map")

    df = po_table()
    isempty(scan.orbits) || append!(df, scan.orbits; cols=:union)
    scan_energies = [m.E for m in scan.maps]

    keys = sort(unique([(o.prime, o.kind) for o in eachrow(df)]))
    section_data = Vector{Dict{Tuple{Int,Int},Vector{Point2f}}}(undef, length(scan_energies))
    boundaries = Vector{Vector{Point2f}}(undef, length(scan_energies))
    background_data = background ? [m.background_points for m in scan.maps] : Vector{Point2f}[]
    info = Vector{String}(undef, length(scan_energies))

    for (i, E) in enumerate(scan_energies)
        section_data[i] = Dict(k => Point2f[] for k in keys)
        at_energy = filter(:E => ==(E), df)
        boundaries[i] = boundary(E; p, n=500)
        info[i] = "E = $(round(E; sigdigits=5)) — $(nrow(at_energy)) periodic orbit(s)"

        for o in eachrow(at_energy)
            key = (o.prime, o.kind)
            append!(section_data[i][key], Point2f.(o.sec_y, o.sec_py))
        end
    end

    all_bounds = reduce(vcat, boundaries; init=Point2f[])
    background && append!(all_bounds, reduce(vcat, background_data; init=Point2f[]))
    if isempty(all_bounds)
        xlim, ylim = (-1.0, 1.0), (-1.0, 1.0)
    else
        xs, ys = getindex.(all_bounds, 1), getindex.(all_bounds, 2)
        xpad = max(0.05 * (maximum(xs) - minimum(xs)), eps())
        ypad = max(0.05 * (maximum(ys) - minimum(ys)), eps())
        xlim = (minimum(xs) - xpad, maximum(xs) + xpad)
        ylim = (minimum(ys) - ypad, maximum(ys) + ypad)
    end

    primes = sort(unique(df.prime))
    kinds = sort(unique(df.kind))

    fig = Figure(size=(1250, 900))
    slider = Slider(fig[2, 1], range=eachindex(scan_energies), startvalue=1)
    idx = slider.value
    ax = Axis(fig[1, 1], xlabel=L"x_2\,[L]", ylabel=L"p_2",
              title=@lift("Periodic orbits, " * info[$idx]))

    if background
        scatter!(ax, @lift(background_data[$idx]); color=(:gray60, 0.5), markersize=2)
    end
    scatter!(ax, @lift(boundaries[$idx]); color=INK, markersize=3)
    for key in keys
        prime, kind = key
        let k = key
            scatter!(ax, @lift(section_data[$idx][k]);
                     color=pick_color(prime), marker=kind_marker(kind), markersize=10)
        end
    end
    xlims!(ax, xlim...)
    ylims!(ax, (-20,20))
    if !isempty(primes) || !isempty(kinds)
        Legend(fig[1, 2],
               [[MarkerElement(color=pick_color(N), marker=:circle, markersize=10) for N in primes],
                [MarkerElement(color=:gray60, marker=kind_marker(k), markersize=10) for k in kinds]],
               [["N = $N" for N in primes], kind_label.(kinds)],
               ["Period", "Stability"];
               framevisible=false, tellheight=false, halign=:left, valign=:top,
               titlehalign=:left, titlesize=12, gridshalign=:left, labelsize=11,
               rowgap=1, groupgap=8)
    end

    on(events(fig.scene).keyboardbutton) do event
        event.action == Keyboard.press || return
        i = slider.value[]
        event.key == Keyboard.right && set_close_to!(slider, min(i + 1, length(scan_energies)))
        event.key == Keyboard.left && set_close_to!(slider, max(i - 1, 1))
    end

    display(fig)
    return (; maps=scan.maps, df, orbits=df, fig, ax, slider, energies=scan_energies)
end
