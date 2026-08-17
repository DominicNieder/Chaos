import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2,
      NonlinearSolve, ADTypes, DataFrames, Dates, ProgressMeter, Printf

include("../styles/makie_theme.jl")
set_theme!(QUARTO_THEME)
include("../models/henon_heiles.jl")
using .HenonHeiles

const CONFIG_DIR = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
const DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
const SAVE_DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/periodic-orbits/")
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
# precompute once, outside the loop
const RGRID  = range(-1.0, 1.0, length = 240)
const EPOT   = [HenonHeiles.potential(x, y, param) for x in RGRID, y in RGRID]
const LEVELS = collect(logrange(5.0*0.009, 6.9*0.089, 7))



"returns -> 1 := elliptic, 2 := hyperbollic (unstable), 3 := hyperbollic (stable), 4 := parabolic,"
function kind_index(τ; ε = 1e-6)
    abs(abs(τ) - 2) < ε && return 4   # parabolic
    abs(τ) < 2 && return 1            # elliptic
    τ > 0 ? 2 : 3                     # hyperbollic (unstable, stable)
end

const KIND_LS    = [:solid, :dash, :dashdot, :dot]
const KIND_MS    = [:circle, :xcross, :diamond, :utriangle]
const KIND_LABEL = ["elliptic", "hyperbolic", "inverse hyperbolic", "parabolic"]


pymax(y, E, p) = sqrt(max(0.0, 2 * p[2] * (E - HenonHeiles.potential(0.0, y, p))))

px2(y, py, E, p)    = 2p[2] * (E - HenonHeiles.potential(0.0, y, p)) - py^2
in_section(v, E, p) = px2(v[1], v[2], E, p) > 0




"Lift (y, py) to a 4D state. Offset follows sign(px) -> no phantom t=0 crossing."
function lift(v, E, p; sgn = +1)
    a = px2(v[1], v[2], E, p)
    a <= 0 && return nothing
    return [EPS_OFF, v[1], sgn * sqrt(a), v[2]]
end

function section_trj(v, prm)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")

    y_sec, py_sec, t_sec = Float64[], Float64[], Float64[]
    cb = HenonHeiles.section_callback(y_sec, py_sec, section_t=t_sec)
    prob = ODEProblem(HenonHeiles.equations!, u0, (0.0, 2000.0), prm.p)
    sol = solve(prob, Vern9(), 
        abstol=1e-14, 
        reltol=1e-14, 
        saveat=dt, 
        callback=cb
        )

    return  sol, permutedims([y_sec py_sec]), t_sec
end

"returns named tuple (;traj, pMap, Nperiod, Tperiod)"
function minPeriodicity(v, prm; tol = 1e-8)
    sol, trace, ts = section_trj(v, prm)
    for i in axes(trace, 2)
        if (norm(v .- trace[:, i]) < tol)   
            k = searchsortedfirst(sol.t, ts[i])
            return (;traj=sol.u[1:k], pMap=trace[:,1:i], Nperiod=i, Tperiod=ts[i])
        end
    end
    return (;traj=sol.u, pMap=trace, Nperiod=nothing, Tperiod=nothing)
end




function creat_integrator(p; tmax = 2000.0, nmax = Ref(1))
    y, py, ts = Float64[], Float64[], Float64[]
    condition(u, t, integ) = u[1]
    affect!(integ) = begin
        push!(y, integ.u[2]); push!(py, integ.u[4]); push!(ts, integ.t)
        length(y) ≥ nmax[] && terminate!(integ)
    end
    cb = ContinuousCallback(condition, affect!, nothing; abstol = 1e-13)
    prob = ODEProblem(HenonHeiles.equations!, zeros(4), (0.0, tmax), p)
    integ = init(prob, Vern9(); abstol = 1e-14, reltol = 1e-14,
                 save_everystep = false, save_start = false, callback = cb)
    return (; integ, y, py, ts, nmax)
end

"prm= (; integ, n, E, p)"
function Tn(v, prm)
    u0 = lift(v, prm.E, prm.p)

    u0 === nothing && error("point $v outside energy boundary")
    # init problem 
    empty!(prm.y); empty!(prm.py); empty!(prm.ts)
    reinit!(prm.integ, u0)
    solve!(prm.integ)

    length(prm.y) < prm.n &&
        error("only $(length(prm.y)) crossings in t < $(prm.tmax) (need $(prm.n))")

    return [prm.y[prm.n], prm.py[prm.n]]
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

function jacobian(v,prm,d)
    J    = zeros(2, 2)
    jacobian!(J,v,prm,d)
end

function find_orbit(v0, prm;
                    N_max = 100, d = 1e-7, tol = 1e-11, 
                    max_backtrack = 30, dmax = 0.05, verbose = false)


    v    = collect(float.(v0))
    hist = zeros(2, N_max)
    J    = zeros(2, 2)

    fail(i, rn, msg) = (v = v, DT = nothing, eigen = nothing, converged = false,
                        history = hist[:, 1:max(i, 0)], resnorm = rn, comment = msg)

    in_section(v, prm.E, prm.p) || return fail(0, Inf, "seed outside energy boundary")

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
            rnew = in_section(vnew, prm.E, prm.p) ? Fres(vnew, prm) : nothing
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

"Is the orbit invariant under S(y,py) = (y,-py)?"
function is_S_symmetric(sec; tol = 1e-6)
    Ssec = sec .* [1.0, -1.0]
    for j in axes(Ssec, 2)
        any(i -> norm(Ssec[:, j] - sec[:, i]) < tol, axes(sec, 2)) || return false
    end
    return true
end


"creating the initial points"
function section_grid(E, p; ny = 10, npy = 10, margin = 0.03)
    roots = real.(filter(r -> abs(imag(r)) < 1e-10, HenonHeiles.limit_of_initial_y0(E, p)))
    sort!(roots)
    ymin, ymax = roots[1], roots[2]      # the two turning points bounding the well
    dy = margin * (ymax - ymin)
    seeds = Vector{Float64}[]
    for y in range(ymin + dy, ymax - dy, length = ny)
        pm = pymax(y, E, p)
        pm <= 0 && continue
        for s in range(0, 1 - margin, length = npy)
            push!(seeds, [y, s * pm])
        end
    end
    seeds
end

function sweep!(df, seeds, prm; tol_min = 1e-8)
    @showprogress dt=1 desc="E=$(round(prm.E,digits=4)) n=$(prm.n)" for v0 in seeds
        res = try
            find_orbit(v0, prm)
        catch e
            @warn "find_orbit threw" seed=v0 exception=e
            continue
        end
        res.converged || continue
        orb = try
            minPeriodicity(res.v, prm; tol = tol_min)
        catch e
            @warn "minPeriodicity threw" v=res.v exception=e
            continue
        end
        orb.Nperiod === nothing && continue

        push!(df, (E = prm.E, n = prm.n,
                   seed_y = v0[1], seed_py = v0[2],
                   y = res.v[1], py = res.v[2],
                   prime = orb.Nperiod, T = orb.Tperiod,
                   trace = tr(res.DT), detDT = det(res.DT),
                   resnorm = res.resnorm, iters = size(res.history, 2),
                   sec = [Vector(c) for c in eachcol(orb.pMap)],
                   history = res.history,
                   traj_x = [u[1] for u in orb.traj],
                   traj_y = [u[2] for u in orb.traj])
                   )
        # add the symmetric orbit, if needed
        sym_py(orb.pMap) && push!(df, (E = prm.E, n = prm.n,
                   seed_y = v0[1], seed_py = v0[2],
                   y = res.v[1], py = -res.v[2],
                   prime = orb.Nperiod, T = orb.Tperiod,
                   trace = tr(res.DT), detDT = det(res.DT),
                   resnorm = res.resnorm, iters = size(res.history, 2),
                   sec = [Vector(c) .* [1.0,-1.0] for c in eachcol(orb.pMap)],
                   history = res.history,
                   traj_x = reverse([u[1] for u in orb.traj]),
                   traj_y = reverse([u[2] for u in orb.traj]))
                   )
    end
    df
end

function run_sweep(Es, ns, p; tmax = 5000.0, ny = 5, npy = 5)
    df = orbit_table()
    timing = DataFrame(E = Float64[], n = Int[], nseeds = Int[],
                       found = Int[], secs = Float64[])
    for E in Es, n in ns
        bff = creat_integrator(p; tmax, nmax = Ref(n))
        prm = (; integ = bff.integ, y = bff.y, py = bff.py, ts = bff.ts,
                 E, p, n, tmax)
        seeds = section_grid(E, p; ny, npy)
        before = nrow(df)
        secs = @elapsed sweep!(df, seeds, prm)
        push!(timing, (E, n, length(seeds), nrow(df) - before, secs))
    end
    (; df, timing)
end


"cleans the data to contain no repetitions"
function dedup(df; tol = 1e-6)
    keep = trues(nrow(df))
    for i in 2:nrow(df)
        for j in 1:i-1
            keep[j] || continue
            if any(norm(a .- b) < tol for a in df.sec[i], b in df.sec[j])
                keep[i] = false
                break
            end
        end
    end
    df[keep, :]
end

"""
    dedup_all(df; tol = 1e-6)

Collapse repeated detections of the same orbit. Orbits are compared only within
the same energy — section points at different E can coincide numerically without
being the same orbit.

Within each energy group, rows are sorted by `prime` so the shortest orbit is
kept as the representative, and its repeats at higher `n` are dropped.
`nhits` records how many raw rows collapsed into each entry (a proxy for basin size).
"""
function dedup_all(df; tol = 1e-6)
    isempty(df) && return copy(df)
    groups = DataFrame[]
    for g in groupby(df, :E)
        push!(groups, dedup(sort(DataFrame(g), [:prime, :resnorm]); tol = tol))
    end
    sort!(vcat(groups...), [:E, :prime, :T])
end

function plot_period(orbits, k; E = 0.1127, bg = nothing, cmap = :managua,
                     show_newton = false, save_fig = false,
                     save_dir = FIG_DIR, show = true)
    sub_E = filter(:E => ==(E), orbits)
    sub_n = filter(:prime => ==(k), sub_E)
    nrow(sub_n) == 0 && return nothing
    cols = resample_cmap(cmap, max(nrow(sub_n), 2))

    # --- configuration space ---
    fc = Figure(size = (1400, 900))
    axc = Axis(fc[1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(),
               title = "config space — prime = $k, $(nrow(sub_n)) orbits, E = $E")
    contour!(axc, RGRID, RGRID, EPOT; levels = LEVELS,
             colormap = :hsv, linewidth = 0.6)
    for (i, row) in enumerate(eachrow(sub_n))
        j = kind_index(row.trace)
        lines!(axc, row.traj_x, row.traj_y;
           color = cols[i], linewidth = 1.2, linestyle = KIND_LS[j])
    end
    Legend(fc[2, 1],
       [LineElement(linestyle = s, color = C_CREAM) for s in KIND_LS],
       KIND_LABEL; orientation = :horizontal, framevisible = false)

    lo, hi = extrema(sub_n.T)
    hi ≈ lo && (hi = lo + max(1e-9, abs(lo)*1e-6))
    Colorbar(fc[1, 2]; limits = (lo, hi), colormap = cmap, label = "T")

    # --- section ---
    fp = Figure(size = (1400, 900))
    axp = Axis(fp[1, 1], xlabel = "y", ylabel = "pᵧ",
               title = "section — prime = $k")
    if bg != nothing
        scatter!(axp, bg...; color = C_RED, markersize = 1.5)
    end
    yb, pb = HenonHeiles.section_boundary_ranges(E, param, 200)
    scatter!(axp, HenonHeiles.section_boundary(yb, pb);
             color = C_CREAM, markersize = 3)
    for (i, row) in enumerate(eachrow(sub_n))
        show_newton && lines!(axp, row.history[1, :], row.history[2, :];
                              color = (C_GOLD, 0.4), linewidth = 0.8)
        j = kind_index(row.trace)
        scatter!(axp, first.(row.sec), last.(row.sec);
                 color = cols[i], markersize = 9, strokewidth = 0.5, marker = KIND_MS[j])
    end
    Legend(fp[2, 1],
       [MarkerElement(marker = m, color =C_CREAM) for m in KIND_MS],
       KIND_LABEL; orientation = :horizontal, framevisible = false)
       
    if save_fig
        mkpath(save_dir)
        tag = tag = "E$(lpad(round(E, digits=4), 6, '0'))_p$(lpad(k, 2, '0'))"
        save(joinpath(save_dir, "config_$tag.png"), fc; px_per_unit = 2)
        save(joinpath(save_dir, "section_$tag.png"), fp; px_per_unit = 2)
    end
    show && (display(GLMakie.Screen(), fc); display(GLMakie.Screen(), fp))
    (; fc, axc, fp, axp, n = nrow(sub_n))

end


function plot_section_all(orbits, E; bg = nothing, ks = nothing,
                          cmap = :managua, label = true,
                          save_fig = false, save_dir = FIG_DIR, show = true)
    sub = filter(row -> row.E ≈ E, orbits)
    nrow(sub) == 0 && return nothing
    ks = something(ks, sort(unique(sub.prime)))
    sub = filter(row -> row.prime in ks, sub)
    cols = resample_cmap(cmap, max(length(ks), 2))
    cidx = Dict(k => i for (i, k) in enumerate(ks))

    fp = Figure(size = (1400, 900))
    axp = Axis(fp[1, 1], xlabel = "y", ylabel = "pᵧ",
               title = "section — E = $(round(E, digits=4)), $(nrow(sub)) orbits")

    bg !== nothing && scatter!(axp, bg...; color = C_RED, markersize = 1.5)
    yb, pb = HenonHeiles.section_boundary_ranges(E, param, 200)
    scatter!(axp, HenonHeiles.section_boundary(yb, pb); color = C_CREAM, markersize = 3)

    for row in eachrow(sub)
        c = cols[cidx[row.prime]]
        j = kind_index(row.trace)
        ys, pys = first.(row.sec), last.(row.sec)
        scatter!(axp, ys, pys; color = c, marker = KIND_MS[j],
                 markersize = 11, strokewidth = 0.5, strokecolor = :black)
        if label
            text!(axp, ys, pys; text = fill(string(row.prime), length(ys)),
                  color = c, fontsize = 11, align = (:left, :bottom),
                  offset = (6, 4))
        end
    end

    Legend(fp[2, 1],
           [MarkerElement(marker = m, color = C_CREAM) for m in KIND_MS],
           KIND_LABEL; orientation = :horizontal, framevisible = false)
    Colorbar(fp[1, 2]; limits = (minimum(ks) - 0.5, maximum(ks) + 0.5),
             colormap = cgrad(cmap, length(ks); categorical = true),
             ticks = ks, label = "prime period")

    if save_fig
        mkpath(save_dir)
        save(joinpath(save_dir, @sprintf("sectionall_E%06.4f.png", E)), fp; px_per_unit = 2)
    end
    show && display(GLMakie.Screen(), fp)
    (; fp, axp, n = nrow(sub))
end

# === Data Types ===

orbit_table() = DataFrame(
    E = Float64[], n = Int[],
    seed_y = Float64[], seed_py = Float64[],
    y = Float64[], py = Float64[],
    prime = Int[], T = Float64[],
    trace = Float64[], detDT = Float64[],
    resnorm = Float64[], iters = Int[],
    sec = Vector{Vector{Float64}}[],
    history = Matrix{Float64}[],
    traj_x = Vector{Float64}[], traj_y=Vector{Float64}[]
)
const Row = @NamedTuple begin
    E::Float64; n::Int
    seed_y::Float64; seed_py::Float64
    y::Float64; py::Float64
    prime::Int; T::Float64
    trace::Float64; detDT::Float64
    resnorm::Float64; iters::Int
    sec::Vector{Vector{Float64}}
    history::Matrix{Float64}
    traj_x::Vector{Float64}; traj_y::Vector{Float64}
end

orbit_table(rows) = DataFrame(filter(!isnothing, rows))


# === here ===
p =(1.0,1.0,1.0)
tmax = 50000.0
max_period=8
E_samples=80
ny, npy = 10, 10

Emin=0.01
Emax=0.1127
Es = collect(range(Emin, Emax, E_samples))
    


println("="^62)
println(" TEST RUN")
println("="^62)

nE_t, np_t, ny_t, npy_t = 2, max_period, 2, 2
test = run_sweep(range(Emin, Emax, nE_t), 1:np_t, p; ny = ny_t, npy = npy_t)

n_full = E_samples * max_period * ny * npy
n_test = nE_t * np_t * ny_t * npy_t
est_min = sum(test.timing.secs) / n_test * n_full / 60

@info "test" raw=nrow(test.df) unique=nrow(dedup_all(test.df)) secs=round(sum(test.timing.secs), digits=1)
@info "estimate" newton_runs=n_full scale="$(round(n_full/n_test, digits=1))×" minutes=round(est_min, digits=1) hours=round(est_min/60, digits=2)

println("="^62)
println("  RUN SWEEP 1")
println("="^62)
out = run_sweep(Es, 1:max_period, p; tmax, ny, npy)
(; df, timing) = out


stamp = Dates.format(now(), "yyyymmdd-HHMM")
jldsave(joinpath(SAVE_DATA_DIR, "orbits_E$(Emin)-$(Emax)_$stamp.jld2");
        df, timing, param, p, tmax, max_period, Es, ny, npy,
        newton_tol = 1e-11, ode_abstol = 1e-14)


orbits = dedup_all(df)

@info "sweep" rows=nrow(df) unique=nrow(orbits) minutes=round(sum(timing.secs)/60, digits=1)
println(combine(groupby(timing, :n), :secs => sum => :secs, :found => sum => :found))
println(combine(groupby(orbits, :prime), nrow => :count))
extrema(orbits.detDT)

# f = joinpath(SAVE_DATA_DIR, "orbits_E0.01-0.1127_20260816-0742.jld2")   
# df,  Es = load(f, "df", "Es")

# orbits = dedup_all(df)           # derived, not saved
# max_period = maximum(df.n)

# for E in Es

#     plot_section_all(orbits, E; bg = nothing, ks = nothing,
#                             cmap = :managua, label = true,
#                             save_fig = true, save_dir = joinpath(FIG_DIR,"20260816-0742B/"), show = false)
# end

# for E in Es, k in 1:max_period
#     plot_period(orbits, k; E, save_fig = true, save_dir=joinpath(FIG_DIR,"20260816-0742/"), show = false)
# end