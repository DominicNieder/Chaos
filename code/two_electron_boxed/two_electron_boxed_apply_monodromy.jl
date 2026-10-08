import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2,
      NonlinearSolve, ADTypes, DataFrames, Dates, ProgressMeter, Printf, ColorSchemes
include("../styles/makie_theme.jl")


# include the integrators
include("integrators-and-setup.jl")




const Row = @NamedTuple begin
    E::Float64 
    v::Vector{Float64}
    T::Float64
    str:: String
end
 
""" 
    E::Float64
    v::Vector{Float64}
    T::Float64
    str:: String
"""
orbit_table() = DataFrame(Row[])





"returns -> 1 := elliptic, 2 := hyperbollic (stable), 3 := hyperbollic (unstable), 4 := parabolic,"
function kind_index(τ; ε = 1e-6)
    abs(abs(τ) - 2) < ε && return 4   # parabolic
    abs(τ) < 2 && return 1            # elliptic
    τ > 0 ? 2 : 3                     # hyperbollic (unstable, stable)
end

const KIND_LS    = [:solid, :dash, :dashdot, :dot]
const KIND_MS    = [:circle, :xcross, :diamond, :utriangle]
const KIND_LABEL = ["elliptic", "hyperbolic", "inverse hyperbolic", "parabolic"]


# =======================================================================
#                       Model 2 electrons in box
# =======================================================================

V_int(u,p)= p.C / abs(u[1] - u[2])  # interaction potential

Kin(u,p)= u[3]^2 / (2p.m1) + u[4]^2 / (2p.m2)  # total kinetic energy

energy(u, p) = Kin(u,p) + V_int(u,p)  # total energy


"""
This system has a upper boundary, for which the energy is not accessible. This is when they the two particles are as far away/close from another as possible. At some point the particle cannot decrease in the system. This is a consequence of an inifit interaction potential V_int, and a finit configuration space. 
"""
function sec_energy_boundry(p)
    sgn1 = p.C > 0 ? 1 : 2          # moving away from the section wall
    sgn2 = p.C > 0 ? 2 : 1
    b1, b2 = get_boxes(p)
    bsec= b1[sgn1]
    ext_pos= b2[sgn2]
    return V_int([bsec, ext_pos, 0, 0], p)
end

# ----------------------------------------------------------------------
# all possible init conditions on manningfold of energy()=E
# ----------------------------------------------------------------------

"Box considered, in units of characteristic box lenght [L]

    (dB(1,1), dB(1,2)), (dB(2,1), dB(2,2))

B1(1|------|2)-B2(1|------|2)"
get_boxes(p::NamedTuple) = ((-(p.L1 + p.del/2), -p.del/2), (p.del/2, (p.L2 + p.del/2)))

function in_box(u::Vector{Float64}, p; tol=0.0)
    b1, b2 = get_boxes(p)
    x1 = u[1]
    x2 = u[2]
    return (b1[1] - tol <= x1 <= b1[2] + tol) &&
    (b2[1] - tol <= x2 <= b2[2] + tol)
end

const EPS_OFF = 1e-15           # offset to the surface of section

section_x1(p) = p.C > 0 ? get_boxes(p)[1][1] + EPS_OFF : get_boxes(p)[1][2] - EPS_OFF

# p1^2 on the section, at the actual wall position
p12(x2, p2, E, p) = 2p.m1 * (E - V_int([section_x1(p), x2, NaN, p2], p) - p2^2/(2p.m2))
in_section(v, E, p) = p12(v[1], v[2], E ,p) > 0

possible_x2(p2,E,p) = 2p.m2*p.C*p.L1/(2p.m2*E-p2^2) + section_x1(p)

x2_is_possible(v,E,p) = v[1] <= possible_x2(v[2],E,p)


"""
    v:: surface of section state
    E:: energy
    p:: system parameters
    return u:: state vector
"""
function lift(v, E, p)
    in_section(v, E, p) || error("v=$v not in section")
    a = p12(v[1], v[2], E, p)
    a > 0 || return nothing
    sgn = p.C > 0 ? 1.0 : -1.0          # moving away from the section wall
    u = [section_x1(p), v[1], sgn*sqrt(a), v[2]]
    in_box(u, p) || error("Did not lift correctly: \nP1: x1=$(u[1]), box1=$(get_boxes(p)[1]);\nP2: x2=$(u[2]), box2=$(get_boxes(p)[2])\nEnergy E(u)=$(energy(u,p)) (set to $(E))\nv=$v")
    return u
end



function analytical_bound(E,p; n=1000)
    dB = get_boxes(p)
    xmax = min(2*p.m2*p.C/p.L2/(2*E*p.m2) + dB[2][1], dB[2][2])
    dx = (xmax - dB[1][2] )/n
    x = collect(dB[1][2]:dx:xmax)
    momentum = sqrt.(max.( (E .- p.C ./ (p.L2 .* (x .+ dB[2][1])) ) .* 2 .* p.m2, 0.0))
    upper = Point2f.(x, momentum)
    lower = Point2f.(reverse(x), -reverse(momentum))
    return vcat(upper, lower, upper[1:1])        # closed curve
end



"""
Calculate the boundary of phase space, where p1=0 and x1=dB(1,j)

    return Point2f.(x2, p), for length of n points, the boundary goes round once
"""
function boundary(E; p, n=400)
    dB = get_boxes(p)
    x1 = section_x1(p)
    lo, hi = dB[2]
    if p.C < 0 && E < 0
        hi = min(hi, x1 + p.C/E)      # attractive: bounded separation
    elseif p.C > 0
        E > 0 || return Point2f[]
        lo = max(lo, x1 + p.C/E)      # repulsive: minimum separation
    end
    lo < hi || return Point2f[]
    x2 = range(lo, hi; length=n)
    p2 = [sqrt(max(2p.m2*(E - V_int([x1, x, 0, 0], p)), 0.0)) for x in x2]
    upper = Point2f.(x2, p2)
    lower = Point2f.(reverse(x2), -reverse(p2))
    return vcat(upper, lower, upper[1:1])
end


"""
    return sol.u, sol.t, pts
"""
function get_traj(u0, t;
    p=(;C=1.0,m1=1.0,m2=1.0, del=1e-3, L1=1.0, L2=1.0), cc_tol=CC_TOL, abstol=INT_TOL, reltol=INT_TOL)
    cb, pts = wall_callback(p; cc_tol=cc_tol)
    prob = ODEProblem(eom!, u0, (0.0, t), p)
    sol  = solve(prob, Vern9(); abstol=abstol, reltol=reltol, callback= cb)
    return sol.u, sol.t, pts
end

"""
    flow ϕₜ takes u(0) to u(t)
    p=(;C, m1, m2, L1, L2,del)
"""
function flow(u0, t; p=(;C= 1.0, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= 1e-8), abstol = INT_TOL, reltol = INT_TOL)
    return get_traj(u0, t; p=p, abstol = abstol, reltol = reltol)[end]
end

function monodromy(u0, t; p=(;C= 1.0, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= 1e-8), d=1e-7, abstol = INT_TOL, reltol = INT_TOL)

    M   = zeros(4, 4)
    for j in 1:4
        e = zeros(4); e[j] = d
        M[:, j] = (flow(u0 .+ e, t; p=p, abstol=abstol, reltol=reltol) .-
                   flow(u0 .- e, t; p=p, abstol=abstol, reltol=reltol)) ./ (2d)
    end
    return M
end

function section_trj(v)
    prm = SectionParams(E, p; 
          tmax = 1_000.0, nfast = 1, ndense=40, 
          cc_tol = CC_TOL, int_tol = INT_TOL, save_everystep = false, save_start = false)
    return section_trj(v, prm)
end

"return sol, pMap, tsd"
function section_trj(v, prm)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")

    empty!(prm.ptsd)
    reinit!(prm.integ_dense, u0)
    solve!(prm.integ_dense)
    sol = prm.integ_dense.sol
    section_points = get_section_p(prm.ptsd)
    trace = isempty(section_points) ? zeros(2, 0) : hcat(section_points...)
    return sol, trace, copy(last.(prm.ptsd))
end

"Find the smallest k with |T^k v - v| < tol.

returns named tuple (;traj, pMap, Nperiod, Tperiod)"
function minPeriodicity(v, prm; pmap_prime_tol = PMAP_PRIME_TOL, search = 40)
    prm.nmax_dense[] = search          # dense integrator only — no interference with T

    sol, trace, ts = section_trj(v, prm)

    for i in axes(trace, 2)
        if norm(v .- trace[:, i]) < pmap_prime_tol
            k = searchsortedfirst(sol.t, ts[i])
            return (; traj = sol.u[1:k], pMap = trace[:, 1:i],
                      Nperiod = i, Tperiod = ts[i])
        end
    end

    @warn "no closure within $search crossings" v pmap_prime_tol
    return (; traj = sol.u, pMap = trace, Nperiod = nothing, Tperiod = nothing)
end





"n-th crossing of the section, starting from v. This is T^n(v)."
function T(v, n::Int, prm::SectionParams)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")
 
    n > prm.nmax_fast[] && (prm.nmax_fast[] = n)
 
    empty!(prm.ptsf)
    reinit!(prm.integ_fast, u0)
    solve!(prm.integ_fast)
 
    length(prm.ptsf) < n &&
        error("only $(length(prm.ptsf)) crossings in t < $(prm.tmax) (need $n)")
    return get_section_p(prm.ptsf)[n]
end

Fres(v, n, prm) = T(v, n, prm) - v


"Finite everywhere: TrustRegion probes outside the boundary and NaN poisons it."
function Fres_safe(v, n, prm::SectionParams)
    a = p12(v[1], v[2], prm.E, prm.p)
    a <= 0 && return fill(1.0 + 100 * sqrt(-a), 2)
    return Fres(v, n, prm)
end


function jacobian!(J, v, n, prm::SectionParams, d)
    E, p = prm.E, prm.p
    for (j, e) in enumerate(([d, 0.0], [0.0, d]))
        vp, vm = v .+ e, v .- e
        okp, okm = in_section(vp, E, p), in_section(vm, E, p)
        if okp && okm
            J[:, j] = (Fres(vp, n, prm) .- Fres(vm, n, prm)) ./ (2d)
        elseif okp
            J[:, j] = (Fres(vp, n, prm) .- Fres(v, n, prm)) ./ d       # forward
        elseif okm
            J[:, j] = (Fres(v, n, prm) .- Fres(vm, n, prm)) ./ d      # backward
        else
            error("both probes outside boundary at $v -- v is on the edge")
        end
    end
    return J
end

function jacobian(v,n,prm,d)
    J    = zeros(2, 2)
    jacobian!(J,v,n,prm,d)
end

function get_DT(v,n,prm; d=1e-7)
    I+jacobian(v,n,prm, d) 
end




function find_orbit(v0, n, prm;
                    N_max = 100, d = 1e-7, tol = PMAP_ROOT_TOL, 
                    max_backtrack = 30, dmax = 0.05, verbose = false)
 
 
    v    = collect(float.(v0))
    hist = zeros(2, N_max)
    J    = zeros(2, 2)
 
    fail(i, rn, msg) = (v = v, DT = nothing, converged = false,
                        history = hist[:, 1:max(i, 0)], resnorm = rn, comment = msg)
 
    in_section(v, prm.E, prm.p) || return fail(0, Inf, "seed outside energy boundary")
 
    r  = Fres(v, n, prm)
    rn = norm(r)
 
    for i in 1:N_max
        hist[:, i] = v
 
        if rn < tol
            jacobian!(J, v, n, prm, d)          # evaluated AT the root
            DT = J + I
            return (v = v, converged = true,
                    history = hist[:, 1:i], resnorm = rn,
                    comment = "|r| = $rn  det(DT) = $(det(DT))")
        end
 
        jacobian!(J, v, n, prm, d)
        κ = cond(J)
        κ > 1e12 && @warn "ill-conditioned Jacobian" iterate=i cond=κ
 
        step = J \ r
        sn   = norm(step)
        sn > dmax && (step .*= dmax / sn)     # cap: keeps Newton local
 
        λ, accepted = 1.0, false              # MUST start false
        for _ in 1:max_backtrack
            vnew = v .- λ .* step
            rnew = in_section(vnew, prm.E, prm.p) ? Fres(vnew, n, prm) : nothing
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


"""
    return (v = sol.u, DT = nothing, converged = true, resnorm = norm(sol.resid),
            history = zeros(2, 0), comment = (sol.retcode)  det(DT) = (det(DT)))

NonlinearSolve variant. AutoFiniteDiff is mandatory: the ODE callback rejects Duals.
"""
function solve_orbit(v0, n, prm; tol = PMAP_ROOT_TOL, maxiters = 300)
    v = collect(float.(v0))
    in_section(v, prm.E, prm.p) ||
        return (v = v, converged = false, resnorm = Inf,
                history = zeros(2, 0), comment = "seed outside boundary")
 
    prob = NonlinearProblem((w, q) -> Fres_safe(w, n, q), v, prm)
    sol  = solve(prob, TrustRegion(autodiff = AutoFiniteDiff(fdtype = Val(:central)));
                 abstol = tol, maxiters)
 
    ok = SciMLBase.successful_retcode(sol) && in_section(sol.u, prm.E, prm.p)
    ok || return (v = sol.u, converged = false,
                  resnorm = norm(sol.resid), history = zeros(2, 0),
                  comment = "$(sol.retcode)")
 
    # DT = get_DT(sol.u, n, prm)
    return (v = sol.u, converged = true, resnorm = norm(sol.resid),
            history = zeros(2, 0), comment = "$(sol.retcode)")
end

"""
    returns a row to the orbit
    return (; E = prm.E, v = res.v, str, T)
"""
function analyse_seed(v0, n, prm; str = "seed",
                      pmap_root_tol = PMAP_ROOT_TOL, pmap_prime_tol= PMAP_PRIME_TOL, maxiters=300)
    # res = find_orbit(v0, n, prm;N_max = 100, d = 1e-7, tol = PMAP_ROOT_TOL, max_backtrack = 30, dmax = 0.05, verbose = false)
    res = solve_orbit(v0, n, prm, tol=pmap_root_tol, maxiters=maxiters)
    # println(res.resnorm, res.comment)
    res.converged || return nothing
    
    min_period = minPeriodicity(res.v, prm; pmap_prime_tol = pmap_prime_tol, search = 40)
    min_period.Nperiod === nothing && return nothing
    T = min_period.Tperiod
    return (;E = prm.E, v = res.v, str, T)
end


function get_eigenvals(M)
    c = tr(M)/2 - 1
    δ = sqrt(Complex(c^2 - 1))
    (c + δ, c - δ)
end


"True if v coincides with ANY section point of an already-catalogued orbit."
function already_found(df, v, prime; pmap_prime_tol = PMAP_PRIME_TOL)
    v === nothing && return false
    for o in eachrow(df)
        o.prime == prime || continue
        any(i -> norm([o.sec_y[i], o.sec_py[i]] .- v) < pmap_prime_tol, eachindex(o.sec_y)) &&
            return true
    end
    return false
end



"creating the initial points from where the search starts from"
function section_grid(E, p; ny = 10, npy = 10, margin = 0.03)
    roots = real.(filter(r -> abs(imag(r)) < 1e-10, HenonHeiles.limit_of_initial_y0(E, p)))
    sort!(roots)
    ymin, ymax = roots[1], roots[2]      # the two turning points bounding the well
    dy = margin * (ymax - ymin)
    seeds = Vector{Float64}[]
    for y in range(ymin + dy, ymax - dy, length = ny)
        pm = pymax(y, E, p)
        pm <= 0 && continue
        for s in range(0 + margin, 1 - margin, length = npy)
            push!(seeds, [y, s * pm])
        end
    end
    seeds
end

"""
Randomly sample the accessible section in `(x2, p2)`.
Positions are sampled uniformly over the accessible `x2` interval, and each
momentum is sampled uniformly over its energy-allowed interval at that `x2`.
"""
function sample_grid(E, p; nx = 10, np = 10, margin = 1e-11)
    nx > 0 || throw(ArgumentError("nx must be positive"))
    np > 0 || throw(ArgumentError("np must be positive"))
    margin >= 0 || throw(ArgumentError("margin must be nonnegative"))

    box2 = get_boxes(p)[2]
    x1 = section_x1(p)
    xmin = box2[1] + margin
    xmax = box2[2] - margin

    # The accessible position interval follows from E - V_int(x1, x2) > 0
    # when p2 = 0. Both boxes lie on opposite sides of the gap, so x2 > x1.
    if p.C > 0
        E > 0 || throw(ArgumentError("no accessible section for C > 0 and E <= 0"))
        xmin = max(xmin, x1 + p.C / E + margin)
    elseif p.C < 0 && E < 0
        xmax = min(xmax, x1 + abs(p.C) / abs(E) - margin)
    elseif p.C == 0
        E > 0 || throw(ArgumentError("no accessible section for C = 0 and E <= 0"))
    end

    xmin < xmax || throw(ArgumentError("the requested energy and margin leave no accessible section"))

    seeds = Vector{Float64}[]
    sizehint!(seeds, nx * np)
    for x2 in xmin .+ rand(nx) .* (xmax - xmin)
        V = V_int([x1, x2, 0.0, 0.0], p)
        p2max_sq = 2p.m2 * (E - V)
        p2max_sq > 0 || continue
        p2max = sqrt(p2max_sq)
        for p2 in (2 .* rand(np) .- 1) .* p2max
            push!(seeds, [x2, p2])
        end
    end
    return seeds
end

"""
    section_map(E, p; nx=5, np=5, tint=1_000, plims=nothing)
 
Seeds the Poincaré section at energy `E` (model units), integrates each seed
for `tint`, and returns:
  - `seed.f`  : seeds on the section with the energy boundary
  - `sec.f`   : surface of section (left), energy drift (right top), legend (right bottom)
  - raw data  : `seeds` (kept), `sec_map`, `trajs`, `times`
"""
function section_map(E, p; nx=5, np=5, tint=1_000, plims=nothing)
    seeds     = sample_grid(E, p; nx, np)
    sec_bound = boundary(E; p, n=500)
 
    # ── integrate ──
    kept    = Vector{Float64}[]
    sec_map = Vector{Point2f}[]
    trajs   = Vector{Vector{Float64}}[]
    times   = Vector{Float64}[]
    for seed in seeds
        u0 = lift(seed, E, p)
        u0 === nothing && continue
        u, ts, pts = get_traj(u0, tint; p)
        push!(kept, seed)
        push!(sec_map, Point2f.(first.(pts, 2)))
        push!(trajs, u)
        push!(times, ts)
    end
 
    # ── seed figure ──
    fseed    = Figure(size=(1300, 900))
    ax_seeds = Axis(fseed[1, 1], xlabel=L"x_2\,[L]", ylabel=L"p_2",
                    title="Seeds, E = $(round(E; sigdigits=4))")
    scatter!(ax_seeds, sec_bound, color=INK, markersize=3)
    scatter!(ax_seeds, Point2f.(kept), color=[pick_color(i) for i in eachindex(kept)],
             markersize=10)
    ylims!(ax_seeds, -20,20)  # if condtion wanted: insert -> plims === nothing || 
 
    # ── section (left) + energy drift and legend (right) ──
    f    = Figure(size=(1300, 900))
    ax   = Axis(f[1:2, 1:2], xlabel=L"x_2\,[L]", ylabel=L"p_2",
                title="Surface of section, E = $(round(E; sigdigits=4))")
    axen = Axis(f[1, 3], xlabel=L"t", ylabel=L"|E(t)-E|\,/\,|E|",
                yscale=log10, title="Energy drift")
 
    Escale = E == 0 ? abs(p.C / p.L1) : abs(E)       # avoid dividing by 0 at E = 0
    scatter!(ax, sec_bound, color=INK, markersize=3)
    for (i, sec) in enumerate(sec_map)
        scatter!(ax, sec, color=pick_color(i), markersize=4)
        err = [max(abs(energy(ui, p) - E) / Escale, eps()) for ui in trajs[i]]
        lines!(axen, times[i], err, color=pick_color(i))
    end
    plims === nothing || ylims!(ax, plims)  #  
 
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
              seeds=kept, sec_map, trajs, times)
end
 
 
"""
    scan_section_maps(Es, p; nx=5, np=5, tint=1_000, plims=nothing,
                      save_fig=false, folder=joinpath(FIG_DIR, "section_maps"))
 
Runs `section_map` for every energy in `Es` (model units). Energies without an
accessible section are skipped with a warning. Returns a vector of results;
view with `display(res[i].sec.f)` or `display(res[i].seed.f)`.
"""
function scan_section_maps(Es, p; nx=5, np=5, tint=1_000, plims=nothing,
                           save_fig=false,
                           folder=joinpath(FIG_DIR, "section_maps_attractive/"))
    save_fig && mkpath(folder)
    tag = @sprintf("C%+g_del%.0e_L%g-%g", p.C, p.del, p.L1, p.L2)
 
    out = []
    @showprogress desc="section maps" for E in Es
        res = try
            section_map(E, p; nx, np, tint, plims)
        catch err
            @warn "skipped E = $E" exception=err
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
#                       Periodic-orbit search up to period nmax
# ============================================================================================

# """
#     E, N (prime period), T (time period), v (root on section),
#     pmap (2×N section points of the orbit), trDT (trace of DT^N), kind (kind_index)
# """
# po_table() = DataFrame(E=Float64[], N=Int[], T=Float64[], v=Vector{Float64}[],
#                        pmap=Matrix{Float64}[], trDT=Float64[], kind=Int[])

# "True if `v` coincides with any section point of an orbit already in `df` at energy `E`."
# function is_known(df, v, E; tol=1e-7)
#     for o in eachrow(df)
#         isapprox(o.E, E; rtol=1e-12, atol=1e-12) || continue
#         any(j -> norm(o.pmap[:, j] .- v) < tol, axes(o.pmap, 2)) && return true
#     end
#     return false
# end


# """
#     find_periodic_orbits(E, p; nmax=8, nx=5, np=5, ...)

# For n = 1…nmax, roots of Tⁿ(v) − v are searched from the same random seeds.
# Each root is reduced to its prime period N (N divides n), deduplicated against
# all section points of already-found orbits, and classified by tr(DTᴺ).
# """
# function find_periodic_orbits(E, p; nmax=8, nx=5, np=5, tmax=100_000.0,
#                               root_tol=PMAP_ROOT_TOL, prime_tol=PMAP_PRIME_TOL,
#                               dup_tol=1e-7, maxiters=300, verbose=false)
#     df    = po_table()
#     seeds = sample_grid(E, p; nx, np)
#     prm   = SectionParams(E, p; tmax, nfast=nmax, ndense=2nmax,
#                           save_everystep=false, save_start=false)

#     for n in 1:nmax
#         prm.nmax_fast[] = n                  # fast integrator stops after n crossings
#         @showprogress desc=@sprintf("E=%.4g  n=%d ", E, n) for v0 in seeds
#             try
#                 res = solve_orbit(v0, n, prm; tol=root_tol, maxiters)
#                 res.converged || continue
#                 is_known(df, res.v, E; tol=dup_tol) && continue

#                 mp = minPeriodicity(res.v, prm; pmap_prime_tol=prime_tol, search=2nmax)
#                 mp.Nperiod === nothing && continue

#                 τ = tr(get_DT(res.v, mp.Nperiod, prm))
#                 push!(df, (E, mp.Nperiod, mp.Tperiod, res.v, mp.pMap, τ, kind_index(τ)))
#             catch err
#                 verbose && @warn "seed failed" E n seed=v0 exception=err
#             end
#         end
#     end
#     return sort!(df, [:N, :T])
# end


# """
#     periodic_orbit_map(E, p; nmax=8, nx=5, np=5, plims=nothing, kwargs...)

# Section with all periodic orbits (colour = period N, marker = stability),
# plus tr(DTᴺ) vs. T on the side (stable band |tr| < 2 shaded).
# """
# function periodic_orbit_map(E, p; nmax=8, nx=5, np=5, plims=nothing, kwargs...)
#     df = find_periodic_orbits(E, p; nmax, nx, np, kwargs...)

#     f    = Figure(size=(1300, 900))
#     ax   = Axis(f[1:2, 1:2], xlabel=L"x_2\,[L]", ylabel=L"p_2",
#                 title="Periodic orbits (N ≤ $nmax), E = $(round(E; sigdigits=4))")
#     axtr = Axis(f[2, 3], xlabel=L"T", ylabel=L"\mathrm{tr}\,DT^N",
#                 yscale=Makie.pseudolog10, title="Stability")

#     scatter!(ax, boundary(E; p, n=500), color=INK, markersize=3)
#     hspan!(axtr, -2, 2, color=(:gray70, 0.3))
#     for o in eachrow(df)
#         scatter!(ax, o.pmap[1, :], o.pmap[2, :], color=pick_color(o.N),
#                  marker=KIND_MS[o.kind], markersize=10)
#         scatter!(axtr, [o.T], [o.trDT], color=pick_color(o.N),
#                  marker=KIND_MS[o.kind], markersize=10)
#     end
#     plims === nothing || ylims!(ax, plims...)

#     if isempty(df)
#         text!(ax, 0.5, 0.5; text="no periodic orbits found", space=:relative,
#               align=(:center, :center))
#     else
#         Ns = sort(unique(df.N)); ks = sort(unique(df.kind))
#         Legend(f[1, 3],
#                [[MarkerElement(color=pick_color(N), marker=:circle, markersize=10) for N in Ns],
#                 [MarkerElement(color=:gray60, marker=KIND_MS[k], markersize=10) for k in ks]],
#                [["N = $N  ($(count(==(N), df.N)))" for N in Ns], KIND_LABEL[ks]],
#                ["Period (count)", "Stability"];
#                framevisible=false, tellheight=false, halign=:left, valign=:top,
#                titlehalign=:left, titlesize=12, gridshalign=:left, labelsize=11,
#                rowgap=1, groupgap=8)
#     end
#     colsize!(f.layout, 3, Auto(0.45))

#     return (; E, df, f, ax, axtr)
# end


# """
#     scan_periodic_orbits(Es, p; nmax=8, nx=5, np=5, plims=nothing,
#                          save_fig=false, save_data=false, kwargs...)

# Runs `periodic_orbit_map` for every energy in `Es` (model units).
# Returns `(; maps, orbits)`: `maps[i].f` to display, `orbits` = all energies in one DataFrame.
# """
# function scan_periodic_orbits(Es, p; nmax=8, nx=5, np=5, plims=nothing,
#                               save_fig=false, save_data=false,
#                               folder=joinpath(FIG_DIR, "periodic_orbits"),
#                               data_dir=DATA_DIR, kwargs...)
#     save_fig  && mkpath(folder)
#     save_data && mkpath(data_dir)
#     tag = @sprintf("C%+g_del%.0e_L%g-%g_N%d", p.C, p.del, p.L1, p.L2, nmax)

#     maps = []
#     for E in Es
#         res = try
#             periodic_orbit_map(E, p; nmax, nx, np, plims, kwargs...)
#         catch err
#             @warn "skipped E = $E" exception=err
#             continue
#         end
#         push!(maps, res)
#         save_fig && save(joinpath(folder, @sprintf("po_E%+.4g_%s.png", E, tag)), res.f)
#     end

#     orbits = isempty(maps) ? po_table() : reduce(vcat, [m.df for m in maps])
#     save_data && jldsave(joinpath(data_dir, "periodic_orbits_$(tag).jld2"); orbits, Es, p)
#     return (; maps, orbits)
# end


include("section_maps.jl")

"""
    Takes a bunch of seed and finds the periodic orbit to each seed. If the orbit is found, it is added to the df. The df is returned at the end.
"""
function sweep!(df::AbstractDataFrame, seeds, n, prm; pmap_root_tol = PMAP_ROOT_TOL, verbose = false)
    @showprogress dt = 1 desc = "E=$(round(prm.E; digits = 4)) n=$n" for v0 in seeds
        orb = try
            analyse_seed(v0, n, prm; str="seed", pmap_root_tol = pmap_root_tol)
        catch e
            verbose && @warn "seed failed" seed = v0 exception = e
            continue
        end
        orb===nothing || push!(df, orb)
    end
    return df
end

function sweep_po_energies(; Es, p=(;C= 1.0, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= 1e-8), nfast = 1, ndense = 2,
    tmax = 100_000.0, ny = 10, np = 10, margin = 0.03, pmap_root_tol = PMAP_ROOT_TOL)
    df = orbit_table()
    for E in Es
        E_inunits = E/ abs(p.C/p.L1)
        prm = SectionParams(E_inunits, p; tmax, nfast, ndense, save_everystep = false, save_start = false)
        seeds = sample_grid(E_inunits, p; nx = ny, np = np, margin = margin)
        sweep!(df, seeds, nfast, prm; pmap_root_tol = pmap_root_tol)
    end
    return df
end

 

"""
finding three base orbits at energy
    return
    (;orb=(; df, timing, seeds), fig=(f, ax))
"""
function get_orbits_1period(;
    p            = (;C= 1.0, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= 1e-8),
    Emin         = 500.0,
    nfast        = 1,              # crossings the dense integrator may take
    ndense       = 2,
    tmax         = 100_000.0,
    seeds = [[0.000154, 20.0], [6.0e-5, 60.0]],  # seeds for the orbits
    orbit_str = ["s", "u"],
    cmap      = COLOR_SCHEME,
    display_figure=false,
    show_figure=true
    )
    Etrue = Emin / abs(p.C/p.L1)
    seeds = sample_grid(Etrue, p; nx = 5, np = 5)
    any(!,[in_section(vi, Etrue, p) for vi in seeds]) && error("some seeds are outside the energy boundary")

    
    prm = SectionParams(Etrue, p; tmax, nfast, ndense, save_everystep = false, save_start = false)

  

    df  = orbit_table()
    timing = DataFrame(E = Float64[], n = Int[], nseeds = Int[],
                       found = Int[], secs = Float64[])

    before = nrow(df)
    secs = @elapsed for i in axes(seeds, 1)
        v0 = seeds[i]
        println("v0=$v0")
        sol = analyse_seed(v0, nfast, prm;str=orbit_str[i])
        sol !== nothing && push!(df, sol)
    end 
    push!(timing, (Etrue, nfast, length(seeds), nrow(df) - before, secs))

    if display_figure  # surface section of the 
        der_rand =  analytical_bound(Etrue,p; n=1000)
        f = Figure(size=(1400,900))
        ax = Axis(f[1, 1], xlabel = "y", ylabel = "py", title = "E=$(round(Etrue, digits=4))")
        vs = Point2f.([o.v[1] for o in eachrow(df)], [o.v[2] for o in eachrow(df)])
        orbit_colors = [pick_color(findfirst(==(o.str), orbit_str), cmap) for o in eachrow(df)]
        scatter!(ax, vs, color = orbit_colors, markersize = 6, alpha = 1)
        annotation = ["$(o.str), T=$(round(o.T, digits=3))" for o in eachrow(df)]
        text!(ax, vs, text=annotation, fontsize = 14, color = INK, align = (:left, :bottom), offset = (5, 5))

        scatter!(ax, Point2f.(vs), markersize = 5, color = INK, alpha = 0.7)
        scatter!(ax, der_rand, markersize = 3, color = INK)
        ylims!(ax, -100,100)
        show_figure && display(f)
        return (;orb=(; df, timing, seeds), fig=(f, ax))
    else
        return (;orb=(; df, timing, seeds), fig=nothing)
    end
end




"""
    follow the orbits A, B, C from low energies upward. The goal is to trace the Monodromy matrix as a function of energy. 
    Start at Energy E=0.11
    return orbs0
"""
function follow_ABC!(orbs0, Es; 
                 p=(;C= 1.0, m1=1.0, m2=1.0, L=1.0, del= 1e-8), tmax = 100_000.0, nfast = 1, ndense = 1, verbose = false)
    prm = SectionParams(orbs0.E[1], p, :both; tmax, nfast, ndense, save_everystep = false, save_start = false)

    for o in collect(eachrow(orbs0))
        v0          = o.v
        orbit_label = o.str
        df          = orbit_table()
        misses      = 0
        @showprogress dt=1 desc="orbit $orbit_label" for e in Es
            set_energy!(prm, e)
            sol = try
                analyse_seed(v0, nfast, prm; str = orbit_label)
            catch err
                verbose && @warn "step failed" orbit=orbit_label E=e exception=err
                nothing
            end
            if sol === nothing
                misses += 1
                misses ≥ 3 && (verbose && @info "orbit $orbit_label lost near E=$e"; break)
                continue
            end
            misses = 0
            push!(df, sol)
            v0 = sol.v
        end
        append!(orbs0, df)
    end
    return (; orb = orbs0)
end

"""
    return (;M, check)

    This takes care of (NaN or Inf) ∈ M
"""
function monodrome(orbits; p = (;C= 1.0, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= 1e-8), verbose=false)
    M    = Vector{Matrix{Float64}}(undef, nrow(orbits))
    check = falses(nrow(orbits))

    @showprogress dt=1 desc="monodromy" for (i, o) in enumerate(eachrow(orbits))
        Mat = try
            monodromy(lift(o.v, o.E, p), o.T; p=p)
        catch err
            verbose && @warn "monodromy failed" E=o.E str=o.str exception=err
            fill(NaN, 4, 4)
        end
        sane     = !(any(isnan, Mat) || any(isinf, Mat)) && abs(det(Mat) - 1.0) < 0.01
        M[i]     = sane ? Mat : fill(NaN, 4, 4)
        check[i] = sane    
    end

    return (; M, check)
end


"
appends .M and .check to orbits

    orbits.mono.M:: Matrix
    orbits.mono.check:: Bool(sensible matrix)
"
function append_monodrome(orbits::DataFrame)
    orb = copy(orbits)
    mono = monodrome(orbits)
    orb.M     = mono.M
    orb.check = mono.check
    return orb
end


function eigen_tr(M)
    τ = tr(M)
    a1= (τ-2)/2
    return a1 .* (1,-1) .* sqrt(a1^2-1)
end

"
    M -> Monodromy matrix
    λ -> max(eigenvalues(M))
    T -> period of periodic orbit
    returns |α|T:: lyapunof exponent to λ= exp(|α|T)
"
get_lyapunov(M::Matrix, T)   = get_lyapunov(eigvals(M), T)
get_lyapunov(λs::AbstractVector, T) = log(maximum(abs.(λs))) 

get_phase2(M::Matrix)          = get_phase2(eigen_tr(M))
get_phase2(λs)                 = angle(λs[2])

# function get_phase(λs::AbstractVector; trivial_tol = 1e-6)
#     transverse = filter(λ -> abs(λ - 1) > trivial_tol, λs)
#     isempty(transverse) && return 0.0
#     λ = argmax(imag, transverse)
#     return angle(λ)
# end
# get_phase(M::Matrix) = get_phase(eigvals(M))

function ABC_energy_trace(;nup=5000,ndown=5000)
    p            = (;C= 1.0, m1=1.0, m2=1.0, L=1.0, del= 1e-8)
    E_fix        = 0.11  # this is the one I fixed
    E_max        = 1.0
    Es_up        = collect(range(E_fix, E_max, nup))[2:end]
    E_min        = 0.01
    Es_down      = sort(collect(range(E_min, E_fix, ndown))[1:end-1]; rev=true)
    nfast        = 1              # crossings the dense integrator may take
    ndense       = 2
    tmax         = 100_000.0
    seeds = [[0.0,0.23], [0.2, 0.3], [0.3, 0.0]]
    orbit_str = ["A", "B", "C"]

    # starting point to follow A, B and C orbit
    res = get_obrits_ABC(p = p, 
                        Emin = E_fix,
                        nfast = nfast, 
                        ndense = ndense, 
                        tmax = tmax,
                        seeds = seeds,
                        orbit_str = orbit_str, 
                        display_figure = false)

    orb_ABC_up, orb_ABC_down=copy(res.orb.df), copy(res.orb.df)  # copy for clean split of variables 


    follow_ABC!(orb_ABC_down,   Es_down)  # follow to lower energy
    follow_ABC!(orb_ABC_up,     Es_up)    # follow to higher energy

    all_ABC = vcat(orb_ABC_down,orb_ABC_up)
    sort!(all_ABC,[:E, :T])

    return (; all_ABC)
end




function graphs(res::DataFrame; fsize=(1400,900), cmap=COLOR_SCHEME)
    df = sort(copy(res), [:E])
    println("="^72)


    function fig_for(label, df)
        sub        = filter(o -> o.str == label, df)

        Ms, check  = ("M" ∈ names(sub) && "check" ∈ names(sub)) ? (sub.M, sub.check) : monodrome(sub)
        Es         = sub.E[check]
        eigmat     = reduce(hcat, eigvals.(Ms[check]))'   # rows = energy, cols = branch j=1:4
        lyp        = map((λs,T) -> get_lyapunov(collect(λs),T), eachrow(eigmat), sub.T[check])
        traces      = [tr(m) for m in Ms[check]]
        θs         = [get_phase2(collect(λs)) for λs in eachrow(eigmat)] ./2pi

        fλ  = Figure(size = fsize)        
        fΘ  = Figure(size = fsize)
        fLy = Figure(size = fsize)
        fTr = Figure(size = fsize)      

        # Lyapunov exponents of (monodrome) Map 
        axLy = Axis(fLy[1,1], xlabel = "E", ylabel = L"|α|T", title = "orbit $label")  
        # trace of Monodromy
        axTr = Axis(fTr[1,1], ylabel = "tr(M)")
        hlines!(axLy,[1.0], color=INK, linewidth=5, linestyle = :dash)
        lines!(axLy, Es, lyp)
        lines!(axTr, Es, traces)

        # Θ plot -> for theta = p//q, p,q∈N, bifurcation happens
        axΘ = Axis(fΘ[1,1], xlabel = "E", ylabel = "Θ/2π", title = "orbit $label", yticklabelcolor = INK, yaxisposition = :left)
        lines!(axΘ, Es, θs)
        
        # eigenvalues λ
        axλ = Axis(fλ[1,1], xlabel = "E", ylabel = "Re(λ) and Im(λ)", title = "orbit $label")
        for j in axes(eigmat, 2)
            lines!(axλ, Es, real.(eigmat[:, j]); color = pick_color(j, cmap), linewidth = 2, linestyle = :solid)
            lines!(axλ, Es, imag.(eigmat[:, j]); color = pick_color(j, cmap), linewidth = 2, linestyle = :dash)

            # text!(ax, Es[end], real(eigmat[end, j]); text = "λ$j",
            #       color = pick_color(j, cmap), align = (:left, :center), offset = (6, 0), fontsize = 13)
        end

        style_elems  = [LineElement(color = :gray70, linestyle = :solid, linewidth = 2),
                        LineElement(color = :gray70, linestyle = :dash,  linewidth = 2)]
        branch_elems = [LineElement(color = pick_color(j, cmap), linewidth = 2) for j in axes(eigmat, 2)]

        Legend(fλ[1,2],
               [style_elems, branch_elems],
               [["Re(λ)", "Im(λ)"], ["λ$j" for j in axes(eigmat, 2)]],
               ["Component", "Branch"])

        (; fλ, axλ, fLy, axLy, fΘ, axΘ, fTr, axTr)
    end

    fA, fB, fC = fig_for("A",df), fig_for("B",df), fig_for("C",df)
    return (; fA , fB , fC)
end

function mask_positive_x(pts::Vector{Point3f})
    out = Point3f[]
    for pt in pts
        if pt[1] <= 0
            push!(out, pt)
        elseif !isempty(out) && !isnan(out[end][1])
            push!(out, Point3f(NaN, NaN, NaN))   # lines! breaks on NaN
        end
    end
    out
end

function halo_lines!(ax, pts; color, linewidth = 3.5, halo_color = (:black,0.5), halo_width = 0.5)
    lines!(ax, pts, color = halo_color, linewidth = linewidth + halo_width)  # dark cover
    lines!(ax, pts, color = color,      linewidth = linewidth)              # colored core
end
"""
    I want a function that generates the the (p,q)-orbits. It is my goal to see or visualize the torus of my preiodic orbits.
"""
function tori(orbits; labels = unique(orbits.str), p = (;C= 1.0, m1=1.0, m2=1.0, L=1.0, del= 1e-8),
               n_periods = 1, n_unst_perido = 10, n_per_shell = 4,
               shell_radius = 3e-3, n_shells = 2,
               colors = COLOR_SCHEME, fig_size = (1400, 1000),
               azimuth = 0.0, elevation = 0.05, perspectiveness = 0.0, mask=false)
 
    df = sort(copy(orbits), [:E])
 
    function torus_of(label)
        sub = filter(o -> o.str == label, df)
        isempty(sub) && error("no orbit found for label $label")
        o = sub[1, :]                        # reference (E, T, v) for this branch
        E, T, v = o.E, o.T, o.v
 
        u0   = lift(v, E, p)
        u0 === nothing && error("orbit $label: seed lies outside the energy boundary")
        core = get_traj(u0, T; p = p, abstol = INT_TOL, reltol = INT_TOL)
 
        shell = Vector{Vector{Float64}}[]
        for k in 0:(n_per_shell - 1), n in 1:n_shells
            θ   = 2π * k / n_per_shell
            δv  = v .+ n .* shell_radius .* [cos(θ), sin(θ)]
            in_section(δv, E, p) || continue
            u0k = lift(δv, E, p)
            u0k === nothing && continue
            if "B" == label
                push!(shell, get_traj(u0k, n_unst_perido * T; p = p, abstol = INT_TOL, reltol = INT_TOL))
            else
                push!(shell, get_traj(u0k, n_periods * T; p = p, abstol = INT_TOL, reltol = INT_TOL))
            end
        end
 
        return (; label, E, T, core, shell)
    end
 
    results = [torus_of(lbl) for lbl in labels]
 
    fig = Figure(size = fig_size)
    # ax  = Axis3(fig[1, 1], xlabel = L"x", ylabel = L"y", zlabel = L"p_y",
    #             title = "Invariant tori around the periodic orbits")
    ax = Axis3(fig[1, 1], xlabel = L"x", ylabel = L"y", zlabel = L"p_y",
               title = "Invariant tori around the periodic orbits",
               azimuth = azimuth, elevation = elevation,
               perspectiveness = perspectiveness)

    for (i, r) in enumerate(results)
        c = colors[mod1(i, length(colors))]
        for traj in r.shell
            pts = mask ?  mask_positive_x(Point3f.(getindex.(traj,1), getindex.(traj,2), getindex.(traj,4))) : Point3f.(getindex.(traj,1), getindex.(traj,2), getindex.(traj,4))
            lines!(ax, pts; color = (c, 1))
        end
        core_pts = Point3f.(getindex.(r.core, 1), getindex.(r.core, 2), getindex.(r.core, 4))
        lines!(ax, core_pts, color = c, linewidth = 4,
               label = "$(r.label)  (E=$(round(r.E, digits=3)), T=$(round(r.T, digits=3)))")
    end
    axislegend(ax, position = :rt)

    return (; fig, ax, results)
end



# res = ABC_energy_trace(nup=5000, ndown=5000).all_ABC

# orbits = append_monodrome(res)


# set_style!(:print)  # :print
# to = tori(orbits; labels=["A", "B", "C"], n_periods=1, 
#             n_unst_perido=1, n_per_shell=3, n_shells=1, 
#             shell_radius = 6e-3, 
#             azimuth = π, elevation = 0.05, perspectiveness = 0.0,
#             mask=false)
# display(to.fig)


# figs = graphs(orbits)


# display(figs.fA.fλ)
# display(figs.fB.fλ)
# display(figs.fC.fλ)

# display(figs.fA.fLy)
# display(figs.fB.fLy)
# display(figs.fC.fLy)

# display(figs.fA.fΘ)
# display(figs.fB.fΘ)
# display(figs.fC.fΘ)


# xlims!(figs.fB.ax, -0.001, 0.1)
# ylims!(figs.fB.ax, -0.5,3)
# save(joinpath(FIG_DIR, "orbitA/Aeigenvalues-vs-E.png"), figs.fA.fig; px_per_unit = 2)
