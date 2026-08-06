#!/usr/bin/env julia
#
# Hénon–Heiles: locate periodic orbits on the x = 0 surface of section and plot them.
#
# Strategy
#   1. Build a grid of admissible seeds (y, p_y) inside the energy boundary.
#   2. ONE pass over the grid: iterate the section map NMAX times per seed and record
#      the return distance ||Q^n(v) - v|| for every n at once.
#   3. Run Newton only from the local minima of that distance map.
#   4. Filter by minimal period, deduplicate, classify stability, plot.
#
# Note on the map Q used here.  The section state is reinjected with p_x = +sqrt(...)
# while the Poincaré map records crossings with direction = -1 (p_x < 0).  Because the
# Hénon–Heiles Hamiltonian is even in x, the reflection S: (x, p_x) -> (-x, -p_x) is a
# symmetry, so this composite map Q = S ∘ P is well defined and Q^2 is the ordinary
# return map on the p_x > 0 branch.  Consequently n counts crossings of x = 0 of EITHER
# sign: an orbit with Q-period n closes after n crossings, i.e. it has ordinary
# Poincaré period n/2 when n is even.  Fixed points of Q^n are genuine periodic points
# in every case, which is all the search needs.
#
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, JSON3, JLD2, Printf

include("../styles/makie_theme.jl")
include("../models/henon_heiles.jl")
using .HenonHeiles

# ----------------------------------------------------------------------
# paths and configuration
# ----------------------------------------------------------------------
CONFIG_FILE = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
DATA_DIR    = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
DATA_FILE   = joinpath(DATA_DIR, "E0.1127-T10000.0-py0.0-n256.jld2")
FIG_DIR     = joinpath(@__DIR__, "../../figures/henon-heiles/periodic-orbits/")
OUT_DIR     = joinpath(@__DIR__, "../../data/henon-heiles/periodic-orbits/")
mkpath(FIG_DIR)
mkpath(OUT_DIR)

configurations = JSON3.read(read(CONFIG_FILE, String))
cfg            = configurations.explore
param          = [Float64(cfg.a.value), Float64(cfg.m.value), Float64(cfg.w.value)]
dt             = Float64(cfg.dt.value)

const E0        = 0.1127
const NMAX      = 1        # longest period searched, in crossings of x = 0
const DY        = 0.02     # scan grid spacing in y
const DPY       = 0.02     # scan grid spacing in p_y
const NEWTON_TOL = 1e-9    # Newton residual tolerance
const DEDUP_TOL  = 1e-4    # fixed points closer than this belong to the same orbit
const PERIOD_TOL = 1e-6    # tolerance for detecting a shorter minimal period

# ----------------------------------------------------------------------
# background: chaotic sea from the stored scan
# ----------------------------------------------------------------------
y_all, py_all = Float32[], Float32[]
if isfile(DATA_FILE)
    data = load(DATA_FILE, "results")
    for d in data
        append!(y_all,  d.sec_y)
        append!(py_all, d.sec_py)
    end
else
    @warn "background scan not found, plotting without it" DATA_FILE
end

y_range, py_range = HenonHeiles.section_boundary_ranges(E0, param, 400)
boundary          = HenonHeiles.section_boundary(y_range, py_range)

# ----------------------------------------------------------------------
# section-map machinery
# ----------------------------------------------------------------------
const PLANE = (1, 0.0)                       # x = 0

"p_x on the section from the energy constraint; NaN outside the accessible region."
function px_from_E(y, py, E, p)
    arg = 2p[2] * (E - HenonHeiles.potential(0.0, y, p)) - py^2
    arg < 0 ? NaN : sqrt(arg)
end

in_section(v, E, p) = !isnan(px_from_E(v[1], v[2], E, p))

"Full phase-space state belonging to a section point, or `nothing` if inadmissible."
function section_state(v, E, p)
    px = px_from_E(v[1], v[2], E, p)
    isnan(px) && return nothing
    # the small offset keeps the root finder from registering t = 0 as a crossing;
    # it caps the attainable residual at roughly 1e-10, hence NEWTON_TOL = 1e-9
    return [1e-10, v[1], px, v[2]]
end

"Q^n(v): section point after n crossings, plus the elapsed time. `nothing` if v is inadmissible."
function Tn(pmap, v, n, E, p)
    u0 = section_state(v, E, p)
    u0 === nothing && return nothing
    reinit!(pmap, u0)
    t0 = current_time(pmap)
    step!(pmap, n)
    u = current_state(pmap)
    return ([u[2], u[4]], current_time(pmap) - t0)
end

"All n crossings as a 2 x n matrix, plus the crossing times."
function iterate_section(pmap, v, n, E, p)
    u0 = section_state(v, E, p)
    u0 === nothing && return nothing
    reinit!(pmap, u0)
    t0    = current_time(pmap)
    trace = zeros(2, n)
    times = zeros(n)
    for i in 1:n
        step!(pmap, 1)
        u           = current_state(pmap)
        trace[:, i] .= (u[2], u[4])
        times[i]    = current_time(pmap) - t0
    end
    return (trace, times)
end

"Residual whose roots are the periodic points."
function Fres(pmap, v, n, E, p)
    out = Tn(pmap, v, n, E, p)
    out === nothing && return nothing
    return out[1] .- v
end

# one system, one map, reused throughout via reinit!
ds   = CoupledODEs(equations!, [-1e-9, 0.0, 0.1, 0.0], param;
                   diffeq = (alg = Vern9(), abstol = 1e-12, reltol = 1e-12))
pmap = PoincareMap(ds, PLANE; direction = +1,
                   rootkw = (xrtol = 1e-12, xatol = 1e-12))

# ----------------------------------------------------------------------
# Newton iteration with finite-difference Jacobian and backtracking
# ----------------------------------------------------------------------
"Central differences where possible, one-sided near the energy boundary."
function jacobian_fd(pmap, v, r0, n, E, p, d)
    J = zeros(2, 2)
    for j in 1:2
        e     = zeros(2); e[j] = d
        rplus  = Fres(pmap, v .+ e, n, E, p)
        rminus = Fres(pmap, v .- e, n, E, p)
        if rplus !== nothing && rminus !== nothing
            J[:, j] = (rplus .- rminus) ./ (2d)
        elseif rplus !== nothing
            J[:, j] = (rplus .- r0) ./ d
        elseif rminus !== nothing
            J[:, j] = (r0 .- rminus) ./ d
        else
            return nothing
        end
    end
    return J
end

function find_orbit(pmap, v0, n, E, p;
                    N_max = 60, d = 1e-6, tol = NEWTON_TOL, max_backtrack = 10)
    v = collect(float.(v0))
    history = Vector{Vector{Float64}}()

    r = Fres(pmap, v, n, E, p)
    r === nothing && return (v = v, converged = false, history = history,
                             resnorm = Inf, J = nothing, comment = "seed outside boundary")
    rn = norm(r)

    for _ in 1:N_max
        push!(history, copy(v))
        if rn < tol
            J = jacobian_fd(pmap, v, r, n, E, p, d)
            return (v = v, converged = true, history = history,
                    resnorm = rn, J = J, comment = @sprintf("residual %.2e", rn))
        end

        J = jacobian_fd(pmap, v, r, n, E, p, d)
        J === nothing && return (v = v, converged = false, history = history,
                                 resnorm = rn, J = nothing,
                                 comment = "Jacobian probe left the boundary")
        abs(det(J)) < 1e-12 && return (v = v, converged = false, history = history,
                                       resnorm = rn, J = J, comment = "singular Jacobian")

        step = J \ r
        # damped step: shrink until the point stays admissible and the residual drops
        accepted = false
        λ = 1.0
        for _ in 1:max_backtrack
            vnew = v .- λ .* step
            rnew = in_section(vnew, E, p) ? Fres(pmap, vnew, n, E, p) : nothing
            if rnew !== nothing && norm(rnew) < rn
                v, r, rn, accepted = vnew, rnew, norm(rnew), true
                break
            end
            λ /= 2
        end
        accepted || return (v = v, converged = false, history = history,
                            resnorm = rn, J = J, comment = "line search stalled")
    end
    return (v = v, converged = false, history = history,
            resnorm = rn, J = nothing, comment = "$N_max iterations exhausted")
end

"Smallest k dividing n with Q^k(v) ≈ v."
function minimal_period(pmap, v, n, E, p; tol = PERIOD_TOL)
    out = iterate_section(pmap, v, n, E, p)
    out === nothing && return n
    trace, _ = out
    for k in 1:(n - 1)
        n % k == 0 && norm(trace[:, k] .- v) < tol && return k
    end
    return n
end

"Residue / linear stability of the periodic point from the map Jacobian DQ^n = J + I."
function classify(J)
    J === nothing && return ("unknown", NaN)
    DQ = J + I
    τ  = tr(DQ)
    kind = abs(τ) < 2 ? "elliptic" : (abs(τ) > 2 ? "hyperbolic" : "parabolic")
    return (kind, τ)
end

# ----------------------------------------------------------------------
# 1. scan pass: return distances for all n in one sweep
# ----------------------------------------------------------------------
ys  = collect(minimum(y_range):DY:maximum(y_range))
pys = collect(-maximum(py_range):DPY:maximum(py_range))

# admissible seeds: inside the energy boundary of the section
admissible = [in_section([y, py], E0, param) for y in ys, py in pys]
@info "seed grid" ny = length(ys) npy = length(pys) admissible = count(admissible)

# --- the seed grid, plotted before any searching happens ---
with_theme(QUARTO_THEME) do
    grid_in  = [(ys[i], pys[j]) for i in eachindex(ys), j in eachindex(pys) if admissible[i, j]]
    grid_out = [(ys[i], pys[j]) for i in eachindex(ys), j in eachindex(pys) if !admissible[i, j]]

    fig = Figure(size = (1400, 950))
    ax  = Axis(fig[1, 1], xlabel = "y", ylabel = "p_y",
               title = @sprintf("Newton seed grid, E = %.4f: %d of %d points admissible (Δy = %.3f, Δp_y = %.3f)",
                                E0, length(grid_in), length(admissible), DY, DPY))

    isempty(y_all) || scatter!(ax, y_all, py_all, color = (:gray, 0.25), markersize = 1.5)
    lines!(ax, boundary, color = C_CREAM, linewidth = 1.5)
    isempty(grid_out) || scatter!(ax, grid_out, color = (:gray, 0.5),
                                  marker = :xcross, markersize = 4)
    scatter!(ax, grid_in, color = C_TEAL, markersize = 5)

    Legend(fig[1, 2],
           [MarkerElement(color = C_TEAL, marker = :circle, markersize = 8),
            MarkerElement(color = (:gray, 0.7), marker = :xcross, markersize = 8)],
           ["admissible", "outside boundary"], "seeds")

    outfig = joinpath(FIG_DIR, @sprintf("seed-grid-E%.4f.png", E0))
    screen = GLMakie.Screen(visible = false)
    display(screen, fig)
    save(outfig, fig, px_per_unit = 1)
    close(screen)
    @info "wrote $outfig"
end

D     = fill(Inf, length(ys), length(pys), NMAX)
valid = falses(length(ys), length(pys))

@info "scanning grid"
for (i, y) in enumerate(ys), (j, py) in enumerate(pys)
    admissible[i, j] || continue
    v   = [y, py]
    out = iterate_section(pmap, v, NMAX, E0, param)
    out === nothing && continue
    trace, _ = out
    valid[i, j] = true
    for n in 1:NMAX
        D[i, j, n] = norm(trace[:, n] .- v)
    end
end

"Grid cells that are local minima of the return distance for this n."
function local_minima(D, valid, n; threshold)
    ny, npy = size(valid)
    seeds = Vector{Tuple{Float64,Int,Int}}()
    for i in 2:(ny - 1), j in 2:(npy - 1)
        valid[i, j] || continue
        d = D[i, j, n]
        (isfinite(d) && d < threshold) || continue
        best = true
        for di in -1:1, dj in -1:1
            (di == 0 && dj == 0) && continue
            if valid[i + di, j + dj] && D[i + di, j + dj, n] < d
                best = false
                break
            end
        end
        best && push!(seeds, (d, i, j))
    end
    sort!(seeds, by = first)
    return seeds
end

# ----------------------------------------------------------------------
# 2. Newton from each candidate, then filter and deduplicate
# ----------------------------------------------------------------------
orbits = Vector{NamedTuple}()   # one entry per distinct periodic orbit
paths  = Dict{Int,Vector{Matrix{Float64}}}()   # Newton paths, keyed by n, for plotting

"Is v already accounted for by one of the orbits found so far?"
function is_new(v, orbits, tol)
    for o in orbits
        for k in axes(o.points, 2)
            norm(o.points[:, k] .- v) < tol && return false
        end
    end
    return true
end

for n in 1:NMAX
    seeds = local_minima(D, valid, n; threshold = 6 * max(DY, DPY))
    paths[n] = Matrix{Float64}[]
    @info "period $n: $(length(seeds)) candidate seeds"

    for (_, i, j) in seeds
        v0  = [ys[i], pys[j]]
        res = find_orbit(pmap, v0, n, E0, param)
        isempty(res.history) || push!(paths[n], reduce(hcat, res.history))

        res.converged || continue
        minimal_period(pmap, res.v, n, E0, param) == n || continue   # drop lower-period repeats
        is_new(res.v, orbits, DEDUP_TOL) || continue

        out = iterate_section(pmap, res.v, n, E0, param)
        out === nothing && continue
        trace, times = out
        kind, τ = classify(res.J)

        push!(orbits, (n = n, v = res.v, points = trace, period = times[end],
                       resnorm = res.resnorm, kind = kind, trace_DQ = τ))
        @printf("n = %d  y = %+.9f  p_y = %+.9f  T = %.6f  |F| = %.2e  %s (tr = %+.4f)\n",
                n, res.v[1], res.v[2], times[end], res.resnorm, kind, τ)
    end
end

@info "found $(length(orbits)) distinct periodic orbits"

# ----------------------------------------------------------------------
# 3. plots
# ----------------------------------------------------------------------
# colour key shared by both figures: one colour per number of crossings
period_colour     = cgrad(:turbo, NMAX; categorical = true)
PLOT_NEWTON_PATHS = true     # set true for per-period diagnostic figures

with_theme(QUARTO_THEME) do
    ns_found = sort(unique(o.n for o in orbits))

    # (a) one Poincaré section carrying every periodic orbit found
    fig = Figure(size = (1600, 1050))
    ax  = Axis(fig[1, 1], xlabel = "y", ylabel = "p_y",
               title = "Hénon–Heiles section x = 0, E = $E0 " *
                       "(circles = elliptic, crosses = hyperbolic)")

    isempty(y_all) || scatter!(ax, y_all, py_all, color = (:gray, 0.3), markersize = 1.5)
    lines!(ax, boundary, color = C_CREAM, linewidth = 1.5)

    for o in sort(orbits, by = x -> x.n)
        scatter!(ax, o.points[1, :], o.points[2, :],
                 color       = period_colour[o.n],
                 marker      = o.kind == "hyperbolic" ? :xcross : :circle,
                 markersize   = 13,
                 strokewidth  = 0.6,
                 strokecolor  = :black)
    end

    if !isempty(ns_found)
        Legend(fig[1, 2],
               [MarkerElement(color = period_colour[n], marker = :circle, markersize = 13)
                for n in ns_found],
               ["n = $n  ($(count(o -> o.n == n, orbits)))" for n in ns_found],
               "crossings")
    end

    outfig = joinpath(FIG_DIR, @sprintf("section-all-periodic-orbits-E%.4f.png", E0))
    screen = GLMakie.Screen(visible = false)
    display(screen, fig)
    save(outfig, fig, px_per_unit = 1)
    close(screen)
    @info "wrote $outfig"

    # (a2) optional: per-period diagnostics showing the Newton paths
    if PLOT_NEWTON_PATHS
        for n in 1:NMAX
            found = filter(o -> o.n == n, orbits)
            (isempty(found) && isempty(get(paths, n, []))) && continue

            fig = Figure(size = (1600, 1000))
            ax  = Axis(fig[1, 1], xlabel = "y", ylabel = "p_y",
                       title = "period-$n points ($(length(found)) orbits)")

            isempty(y_all) || scatter!(ax, y_all, py_all, color = (:gray, 0.3), markersize = 1.5)
            lines!(ax, boundary, color = C_CREAM, linewidth = 1.5)

            for h in paths[n]
                scatterlines!(ax, h[1, :], h[2, :], color = C_PURPLE,
                              markersize = 4, linewidth = 1)
            end
            for o in found
                scatter!(ax, o.points[1, :], o.points[2, :],
                         marker = :cross, markersize = 18, color = C_TEAL)
            end

            outfig = joinpath(FIG_DIR, @sprintf("section-period%02d.png", n))
            screen = GLMakie.Screen(visible = false)
            display(screen, fig)
            save(outfig, fig, px_per_unit = 1)
            close(screen)
        end
    end

    # (b) all orbits in configuration space, one figure, colour = period
    fig = Figure(size = (1100, 950))
    ax  = Axis(fig[1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(),
               title = "Hénon–Heiles periodic orbits, E = $E0 " *
                       "(solid = elliptic, dashed = hyperbolic)")

    # equipotential V(x, y) = E, the boundary of the accessible region
    xs_c  = range(-1.1, 1.1, length = 500)
    ys_c  = range(-0.8, 1.3, length = 500)
    Vgrid = [HenonHeiles.potential(x, y, param) for x in xs_c, y in ys_c]
    contour!(ax, xs_c, ys_c, Vgrid, levels = [E0],
             color = C_CREAM, linewidth = 1.5)

    for o in sort(orbits, by = x -> x.n)
        u0  = [1e-10, o.v[1], px_from_E(o.v[1], o.v[2], E0, param), o.v[2]]
        dso = CoupledODEs(equations!, u0, param;
                          diffeq = (alg = Vern9(), abstol = 1e-12, reltol = 1e-12))
        X, _ = trajectory(dso, o.period, u0; Δt = min(dt, o.period / 4000))

        lines!(ax, X[:, 1], X[:, 2],
               color     = (period_colour[o.n], 0.85),
               linewidth = 2,
               linestyle = o.kind == "hyperbolic" ? :dash : :solid)
    end

    # legend: one entry per period actually present
    if !isempty(ns_found)
        Legend(fig[1, 2],
               [LineElement(color = period_colour[n], linewidth = 3) for n in ns_found],
               ["n = $n  ($(count(o -> o.n == n, orbits)))" for n in ns_found],
               "crossings")
    end

    outfig = joinpath(FIG_DIR, @sprintf("all-periodic-orbits-E%.4f.png", E0))
    screen = GLMakie.Screen(visible = false)
    display(screen, fig)
    save(outfig, fig, px_per_unit = 1)
    close(screen)
    @info "wrote $outfig"
end

# ----------------------------------------------------------------------
# 4. save
# ----------------------------------------------------------------------
jldsave(joinpath(OUT_DIR, @sprintf("periodic-orbits-E%.4f-n%d.jld2", E0, NMAX));
        orbits = orbits, E = E0, param = param)

open(joinpath(OUT_DIR, @sprintf("periodic-orbits-E%.4f-n%d.json", E0, NMAX)), "w") do io
    JSON3.pretty(io, [Dict("n" => o.n, "y" => o.v[1], "py" => o.v[2],
                           "period" => o.period, "resnorm" => o.resnorm,
                           "stability" => o.kind, "trace" => o.trace_DQ) for o in orbits])
end