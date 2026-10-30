using Roots

uc(c) = find_zero(u -> 2/π*(atan(u) + u/(1+u^2)) - c, (0.0, 1e6), Roots.Bisection())
ε(E, m, c) = uc(c) * sqrt(-2m*E)


"""
    section_map(E, p; nx=5, np=5, tint=1_000, plims=nothing, parallel=true)

Seeds the Poincaré section at energy `E` (model units), integrates each seed
for `tint`, and returns:
  - `seed.f`  : seeds on the section with the energy boundary
  - `sec.f`   : surface of section (left), energy drift (right top), legend (right bottom)
  - raw data  : `seeds` (kept), `sec_map`

Seed trajectories are integrated in parallel when Julia is running with
multiple threads; plotting is always performed serially.
"""
function section_map(E, p; nx=5, np=5, tint=1_000, plims=(-ε(E, 1, 0.99),ε(E, 1, 0.99)), vis=:dark, set_pmin=nothing, set_pmax=nothing, set_xmin=nothing, set_xmax=nothing, parallel=true, show_progress=true)
    set_style!(vis)
    seeds     = sample_grid(E, p; nx, np, set_pmin=set_pmin, set_pmax=set_pmax, set_xmin=set_xmin, set_xmax=set_xmax)
    sec_bound = boundary(E; p, n=500)

    f    = Figure(size=(1300, 900))
    ax   = Axis(f[1:2, 1:2], xlabel=L"x_2\,[L]", ylabel=L"p_2",
                title="Surface of section, E = $(round(E; sigdigits=4))")
    axen = Axis(f[1, 3], xlabel=L"t", ylabel=L"|E(t)-E|\,/\,|E|",
                yscale=log10, title="Energy drift")
    Escale = E == 0 ? abs(p.C / p.L1) : abs(E)
    scatter!(ax, sec_bound, color=INK, markersize=3)

    # ── integrate independently per seed; keep plotting serial ──
    results = Vector{Union{Nothing,NamedTuple}}(nothing, length(seeds))
    function integrate_seed(i)
        seed = seeds[i]
        u0 = lift(seed, E, p)
        u0 === nothing && return nothing
        u, ts, pts = get_traj(u0, tint; p)
        err = [max(abs(energy(ui, p) - E) / Escale, eps()) for ui in u]
        return (; seed, u, ts, pts, err)
    end

    if parallel && Threads.nthreads() > 1 && length(seeds) > 1
        @showprogress enabled=show_progress desc="seed trajectories" Threads.@threads for i in eachindex(seeds)
            results[i] = integrate_seed(i)
        end
    else
        parallel && length(seeds) > 1 && Threads.nthreads() == 1 &&
            @warn "Julia is running with one thread; start Julia with --threads=auto to parallelize seed trajectories"
        @showprogress enabled=show_progress desc="seed trajectories" for i in eachindex(seeds)
            results[i] = integrate_seed(i)
        end
    end

    kept    = Vector{Float64}[]
    sec_map = Vector{Point2f}[]

    naborted = 0
    for result in results
        result === nothing && continue
        seed, u, ts, pts, err = result.seed, result.u, result.ts, result.pts, result.err
        # if ts[end] < tint                      # solver aborted (near-collision, dt → eps)
        #     naborted += 1
        #     continue
        # end
        push!(kept, seed)
        sec = Point2f.(first.(pts, 2))
        push!(sec_map, sec)
        scatter!(ax, sec, color=pick_color(length(kept)), markersize=4)
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
                      save_fig=false, parallel=true, show_progress=true,
                      folder=joinpath(FIG_DIR, "section_maps"))

Runs `section_map` for every energy in `Es` (model units). Energies without an
accessible section are skipped with a warning. Returns a vector of results;
view with `display(res[i].sec.f)` or `display(res[i].seed.f)`. Each map's seed
trajectories are integrated in parallel when multiple Julia threads are
available; map generation and plotting remain serial.
"""
function scan_section_maps(Es, p; nx=5, np=5, tint=1_000,
                           save_fig=false,
                           folder=joinpath(FIG_DIR, "section_maps"), c=0.99,
                           parallel=true, show_progress=true)
    save_fig && mkpath(folder)
    tag = @sprintf("C%+g_del%.0e_L%g-%g", p.C, p.del, p.L1, p.L2)

    out = []
    @showprogress enabled=show_progress desc="section maps" for E in Es
        res = try
            plims = (-ε(E, 1.0, c), ε(E, 1.0, c))
            section_map(E, p; nx, np, tint, plims, parallel, show_progress=false)
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
(stable band |τ| < 2 shaded). Set `parallel=false` to disable threaded
integration of background-map seeds; plotting remains serial.
"""
function periodic_orbit_map(E, p; nmax=8, nx=5, np=5, plims=(-ε(E, 1, 0.99),ε(E, 1, 0.99)),
                            background=true, nbg=3, tint=1_000, orbit_df=nothing,
                            parallel=true, show_progress=true, kwargs...)
    df = orbit_df === nothing ?
         find_periodic_orbits(E, p; nmax, nx, np, show_progress, kwargs...) : orbit_df
    background_points = Point2f[]

    f    = Figure(size=(1300, 900))
    ax   = Axis(f[1:2, 1:2], xlabel=L"x_2\,[L]", ylabel=L"p_2",
                title="Periodic orbits (N ≤ $nmax), E = $(round(E; sigdigits=4))")
    axtr = Axis(f[2, 3], xlabel=L"T", ylabel=L"\tau = \mathrm{tr}\,M - 2",
                yscale=Makie.pseudolog10, title="Stability")

    if background
        bg = section_map(E, p; nx=nbg, np=nbg, tint, parallel, show_progress)
        for sec in bg.sec_map
            scatter!(ax, sec, color=(:gray60, 0.7), markersize=2)
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
are available, then builds the Makie maps serially (including any background
section maps; their seed integrations can also run in parallel). Set
`parallel=false` to disable threaded work. Progress is reported during the
energy scan and background integrations; set `show_progress=false` to suppress
it. Start Julia with `--threads=auto` to use all available CPU cores.
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
                               background, nbg, tint, orbit_df=orbit_dfs[i],
                               parallel, show_progress, kwargs...)
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



# ============================================================================================
#        Interactive section-map scan (no periodic-orbit search)
#        Replaces the empty `interactive_section_maps` stub in section_maps.jl
# ============================================================================================

"""
    relative_seeds(E, p; nx=5, np=5, margin=0.02)

Deterministic seeds at fixed *relative* positions in the accessible section:
x₂ at fraction `sx ∈ [margin, 1-margin]` of the accessible x₂-interval, p₂ at
fraction `sp ∈ [-1+margin, 1-margin]` of p₂max(x₂). Seed k therefore "means the same"
at every energy, so its section curve can be followed continuously through `Es`.
Inaccessible slots are `nothing` (slot count is always nx*np).
"""
function relative_seeds(E, p; nx=5, np=5, margin=0.02)
    seeds = Vector{Union{Nothing,Vector{Float64}}}(nothing, nx*np)
    lo, hi = get_boxes(p)[2]
    x1 = section_x1(p)
    if p.C > 0
        E > 0 || return seeds
        lo = max(lo, x1 + p.C/E)
    elseif p.C < 0 && E < 0
        hi = min(hi, x1 + abs(p.C)/abs(E))
    end
    hi > lo || return seeds
    w   = hi - lo
    sxs = nx == 1 ? [0.5] : range(margin, 1 - margin, length=nx)
    sps = np == 1 ? [0.0] : range(-1 + margin, 1 - margin, length=np)
    k = 0
    for sx in sxs, sp in sps
        k += 1
        x2   = lo + sx*w
        p2sq = 2p.m2*(E - V_int([x1, x2, 0.0, 0.0], p))
        p2sq > 0 || continue
        seeds[k] = [x2, sp*sqrt(p2sq)]
    end
    return seeds
end

"Integrate every seed at energy E; returns per-slot section points and max relative energy drift."
function compute_section(E, p, seeds; tint=1_000)
    secs   = [Point2f[] for _ in seeds]
    drift  = fill(NaN, length(seeds))
    Escale = E == 0 ? abs(p.C/p.L1) : abs(E)
    for (k, s) in enumerate(seeds)
        s === nothing && continue
        try
            in_section(s, E, p) || continue
            u0 = lift(s, E, p)
            u0 === nothing && continue
            u, _, pts = get_traj(u0, tint; p)
            secs[k]  = Point2f.(first.(pts, 2))
            drift[k] = max(maximum(ui -> abs(energy(ui, p) - E), u) / Escale, eps())
        catch err
            @debug "seed failed" E seed=s exception=err
        end
    end
    return secs, drift
end

"""
    interactive_section_maps(Es, p; nx=5, np=5, tint=1_000, seeds=:relative,
                             plims=nothing, follow_limits=true, vis=:dark,
                             markersize=3, parallel=true, show_progress=true,
                             folder=joinpath(FIG_DIR, "section_maps"))

Computes the surface of section for every energy in `Es` (model units) and shows
them in one GLMakie window with an energy slider. No periodic-orbit search.

`seeds`:
  - `:relative` (default) — `relative_seeds`, same relative position at each E,
                            so each colour can be followed through the scan
  - `:random`             — `sample_grid` (re-randomised at each E)
  - a vector of `[x₂, p₂]` — fixed absolute seeds, skipped where inaccessible

Controls: slider, ←/→ step through energies, `s` saves the current frame to `folder`.
Right panel: max relative energy drift vs E (worst seed), with the current E marked.
`follow_limits=true` rescales the axes to each energy; `plims` fixes the p₂ range.

Returns `(; fig, ax, axdrift, slider, energies, secs, drift, seeds)`.
"""
function interactive_section_maps(Es, p; nx=5, np=5, tint=1_000, seeds=:relative,
                                  plims=nothing, follow_limits=true, vis=:dark,
                                  markersize=3, parallel=true, show_progress=true,
                                  folder=joinpath(FIG_DIR, "section_maps"))
    set_style!(vis)
    energies = collect(float.(Es))
    nE = length(energies)

    seed_at(E) = seeds === :relative ? relative_seeds(E, p; nx, np) :
                 seeds === :random   ? (try sample_grid(E, p; nx, np) catch; Vector{Float64}[] end) :
                 [in_section(s, E, p) ? s : nothing for s in seeds]

    seedsets = [seed_at(E) for E in energies]
    nslots   = maximum(length, seedsets; init=0)
    nslots > 0 || error("no accessible seeds at any energy")

    secs  = Vector{Vector{Vector{Point2f}}}(undef, nE)
    drift = Vector{Vector{Float64}}(undef, nE)

    work(i) = ((secs[i], drift[i]) = compute_section(energies[i], p, seedsets[i]; tint))
    if parallel && Threads.nthreads() > 1 && nE > 1
        @showprogress enabled=show_progress desc="section maps" Threads.@threads for i in 1:nE
            work(i)
        end
    else
        parallel && nE > 1 && Threads.nthreads() == 1 &&
            @warn "one Julia thread; start with --threads=auto to parallelize"
        @showprogress enabled=show_progress desc="section maps" for i in 1:nE
            work(i)
        end
    end

    bounds   = [boundary(E; p, n=500) for E in energies]
    maxdrift = [all(isnan, d) ? NaN : maximum(filter(!isnan, d)) for d in drift]
    nkept    = [count(!isempty, s) for s in secs]

    # ── figure ──
    fig    = Figure(size=(1300, 900))
    slider = Slider(fig[2, 1:2], range=1:nE, startvalue=1)
    idx    = slider.value

    ax = Axis(fig[1, 1], xlabel=L"x_2\,[L]", ylabel=L"p_2",
              title=@lift("Surface of section, E = $(round(energies[$idx]; sigdigits=5))  " *
                          "($(nkept[$idx]) / $(length(seedsets[$idx])) seeds,  $(idx[]) / $nE)"))
    axdrift = Axis(fig[1, 2], xlabel=L"E", ylabel=L"\max_t |E(t)-E|/|E|",
                   yscale=log10, title="Energy drift")

    scatter!(ax, @lift(bounds[$idx]); color=INK, markersize=3)
    for k in 1:nslots                                  # one plot per slot: no length mismatch
        scatter!(ax, @lift(k <= length(secs[$idx]) ? secs[$idx][k] : Point2f[]);
                 color=pick_color(k), markersize=markersize)
    end

    good = .!isnan.(maxdrift)
    any(good) && scatterlines!(axdrift, energies[good], maxdrift[good];
                               color=INK, markersize=4)
    vlines!(axdrift, @lift([energies[$idx]]); color=pick_color(1), linewidth=2)
    colsize!(fig.layout, 2, Relative(0.3))

    # ── limits ──
    function set_limits!(i)
        follow_limits || return
        pts = vcat(bounds[i], reduce(vcat, secs[i]; init=Point2f[]))
        pts = filter(q -> all(isfinite, q), pts)
        isempty(pts) && return
        xs, ys = first.(pts), last.(pts)
        xp = max(0.03*(maximum(xs) - minimum(xs)), eps(Float32))
        yp = max(0.03*(maximum(ys) - minimum(ys)), eps(Float32))
        xlims!(ax, minimum(xs) - xp, maximum(xs) + xp)
        plims === nothing ? ylims!(ax, minimum(ys) - yp, maximum(ys) + yp) : ylims!(ax, plims...)
    end
    on(set_limits!, idx)
    set_limits!(1)
    !follow_limits && plims !== nothing && ylims!(ax, plims...)

    # ── keyboard: ←/→ step, s saves the current frame ──
    on(events(fig.scene).keyboardbutton) do ev
        ev.action == Keyboard.press || return
        i = slider.value[]
        ev.key == Keyboard.right && set_close_to!(slider, min(i + 1, nE))
        ev.key == Keyboard.left  && set_close_to!(slider, max(i - 1, 1))
        if ev.key == Keyboard.s
            mkpath(folder)
            tag  = @sprintf("C%+g_del%.0e_L%g-%g", p.C, p.del, p.L1, p.L2)
            file = joinpath(folder, @sprintf("isec_E%+.5g_%s.png", energies[i], tag))
            save(file, fig); @info "saved $file"
        end
    end

    display(fig)
    return (; fig, ax, axdrift, slider, energies, secs, drift, seeds=seedsets)
end



function sec_map_vid500_2(;n=64)
    Es = range(-500,-2, n)
    p= (;C=-1, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= 1e-8)
    seedE0 = [[]]  # select seeds for E=-500, then stretch x by how the ration increased
    for e in Es
        ε(E, p.m2, 0.99)
    end
    interactive_section_maps(Es, p; nx=5, np=5, tint=1_000, seeds=:relative,
                                  plims=nothing, follow_limits=true, vis=:dark,
                                  markersize=3, parallel=true, show_progress=true,
                                  folder=joinpath(FIG_DIR, "section_maps"))
end