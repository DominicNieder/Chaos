import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include("a_tester.jl")


using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2,
    NonlinearSolve, ADTypes, DataFrames, Dates, ProgressMeter, Printf, ColorSchemes,
    Test

include("../styles/makie_theme.jl")
set_theme!(QUARTO_THEME)
include("../models/henon_heiles.jl")
using .HenonHeiles

BLAS.set_num_threads(1)   # all linear algebra here is 2x2; BLAS threads only compete


# =====================================================================
#  Paths and configuration
# =====================================================================

const CONFIG_DIR    = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
const DATA_DIR      = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
const FIG_DIR       = joinpath(@__DIR__, "../../figures/henon-heiles/periodic-orbits/")
const SAVE_DATA_DIR = joinpath(@__DIR__, "../../data/henon-heiles/periodic-orbits/")
const DATA_FILE     = joinpath(DATA_DIR, "E0.1127-T10000.0-py0.0-n256.jld2")




# =====================================================================
#  Tolerances
# =====================================================================

const EPS_OFF       = 1e-9    # offset from the section; must exceed TOL_CALLBACK
const TOL_INTEGRATE = 1e-14   # ODE abstol = reltol
const TOL_CALLBACK  = 1e-13   # root finding on the section; looser than the integrator
const TOL_PO_ROOT   = 1e-11   # Newton convergence, |T^n(v) - v|
const TOL_PERIOD    = 1e-8    # closure test in minPeriodicity
const TOL_ROOT_COMP = 1e-6    # two section points are "the same" if closer than this
const FD_STEP       = 1e-7    # finite-difference step for the Jacobian
const TMIN_CROSS    = 1e-6    # ignore section crossings this soon after reinit


# =====================================================================
#  Plot constants
# =====================================================================
const DT     = 0.3      
const PARAM  = (1.0, 1.0,1.0)
const RGRID  = range(-1.0, 1.0, length = 240)
const EPOT   = [HenonHeiles.potential(x, y, PARAM) for x in RGRID, y in RGRID]
const LEVELS = collect(logrange(5.0 * 0.009, 6.9 * 0.089, 7))

const KIND_LS    = [:solid, :dash, :dashdot, :dot]
const KIND_MS    = [:circle, :xcross, :diamond, :utriangle]
const KIND_LABEL = ["elliptic", "hyperbolic", "inverse hyperbolic", "parabolic"]

const COLOR_SCHEME = [C_ORANGE, C_TEAL, C_CREAM, C_GOLD, C_PURPLE, C_GREEN]
pick_color(i) = COLOR_SCHEME[mod1(i, length(COLOR_SCHEME))]

"1 := elliptic, 2 := hyperbolic, 3 := inverse hyperbolic, 4 := parabolic"
function kind_index(τ; ε = 1e-6)
    abs(abs(τ) - 2) < ε && return 4
    abs(τ) < 2          && return 1
    τ > 0 ? 2 : 3
end


# =====================================================================
#  The symmetry group  G = D3 x Z2,  order 12
# =====================================================================
#
#  V(x, y) = m w^2 (x^2 + y^2)/2 + a (x^2 y - y^3/3)
#          = m w^2 r^2 / 2 + a r^3 sin(3θ) / 3
#
#  so the potential is invariant under
#      C3 : θ -> θ + 2π/3          (rotation)
#      σ  : x -> -x                (reflection in the y-axis; θ -> π - θ)
#  which generate D3, order 6. Note that y -> -y alone is NOT a symmetry:
#  the cubic term is odd in y at fixed x.
#
#  Time reversal Θ: p -> -p, t -> -t, is anti-canonical and commutes with D3,
#  giving G = D3 x Z2 of order 12. Every g in G preserves E, the period T and
#  tr(DT), so a whole G-orbit of periodic orbits carries identical dynamical
#  data: one representative per class is enough for analysis.
#
#  An element is stored as (name, M, τ) with M the 2x2 configuration-space
#  matrix and τ = ±1 the time-reversal flag. It acts on phase space as
#      (q, p)  ->  (M q, τ M p).

const SYM_GROUP = let
    R  = [cospi(2/3) - sinpi(2/3); sinpi(2/3) cospi(2/3)]
    S  = [-1.0 0.0; 0.0 1.0]
    d3 = [("E", Matrix{Float64}(I, 2, 2)), ("C3", R), ("C3²", R * R)]
    d3 = vcat(d3, [(n == "E" ? "σ" : "σ" * n, S * M) for (n, M) in d3])
    [(; name = τ == 1 ? n : n * "Θ", M = M, τ = τ) for τ in (1, -1) for (n, M) in d3]
end

"Action of a group element on a 4D phase-space state."
function apply_sym(g, u)
    q = g.M * [u[1], u[2]]
    p = g.τ .* (g.M * [u[3], u[4]])
    return [q[1], q[2], p[1], p[2]]
end

"Group multiplication, by lookup. Doubles as a closure check."
function compose(g, h; tol = 1e-9)
    M, τ = g.M * h.M, g.τ * h.τ
    for k in SYM_GROUP
        k.τ == τ && norm(k.M - M) < tol && return k
    end
    error("composition $(g.name)*$(h.name) fell outside the group")
end

"""
Verify at load time that every element really is a symmetry: M orthogonal (so
the map is canonical up to the time flip), the potential invariant, and the
multiplication table closed.
"""
function check_group(p = PARAM; n = 200, tol = 1e-11, rng = MersenneTwister(1))
    for g in SYM_GROUP
        norm(g.M' * g.M - I) < 1e-12 || error("$(g.name): M is not orthogonal")
        for _ in 1:n
            x, y = 2rand(rng) - 1, 2rand(rng) - 1
            q = g.M * [x, y]
            abs(HenonHeiles.potential(q[1], q[2], p) -
                HenonHeiles.potential(x, y, p)) < tol ||
                error("$(g.name) does not preserve the potential")
        end
        for h in SYM_GROUP
            compose(g, h)          # throws if the table is not closed
        end
    end
    return true
end

"""
Does g map the section {x = 0, px > 0} to itself?

M must fix the y-axis pointwise and px must come back positive. Only two
elements qualify: the identity, and σΘ, which acts as (y, py) -> (y, -py).
Everything else leaves the section and needs one integration to come back.
"""
function is_section_stabiliser(g; tol = 1e-12)
    norm(g.M * [0.0, 1.0] - [0.0, 1.0]) < tol || return false
    return g.τ * (g.M * [1.0, 0.0])[1] > 0
end

"Action of a section-stabilising element on (y, py). M fixes the y-axis, so only τ acts."
sec_action(g, v) = [v[1], g.τ * v[2]]

const SEC_STAB = filter(is_section_stabiliser, SYM_GROUP)

"One representative per coset g*SEC_STAB. |G|/|SEC_STAB| = 6 flows per orbit."
const COSET_REPS = let
    reps, seen = eltype(SYM_GROUP)[], Set{String}()
    for g in SYM_GROUP
        g.name in seen && continue
        push!(reps, g)
        for h in SEC_STAB
            push!(seen, compose(g, h).name)
        end
    end
    reps
end

check_group()
@assert length(SYM_GROUP)  == 12
@assert length(SEC_STAB)   == 2
@assert length(COSET_REPS) == 6


# =====================================================================
#  Data structures
# =====================================================================

struct SectionParams{IF, ID, P}
    integ_fast  :: IF
    yf          :: Vector{Float64}
    pyf         :: Vector{Float64}
    tsf         :: Vector{Float64}
    nmax_fast   :: Base.RefValue{Int}

    integ_dense :: ID
    yd          :: Vector{Float64}
    pyd         :: Vector{Float64}
    tsd         :: Vector{Float64}
    nmax_dense  :: Base.RefValue{Int}

    E           :: Float64
    p           :: P
    tmax        :: Float64
end

"""
One periodic orbit.

`class`/`index` are the linear stability class; the `sym_*` fields describe the
symmetry class. `sym_y`/`sym_py` are section points of the images of this orbit
under all 12 group elements, `sym_stab` lists the elements that map the orbit to
itself, and `sym_mult` = 12 / |stabiliser| is the number of distinct orbits in
the class. `sym_id` is assigned by `dedup_classes` and is shared by every row of
one class at one energy.

The symmetry fields stay empty until `annotate_symmetry!` has run.
"""
const Row = @NamedTuple begin
    E::Float64; n::Int
    seed_y::Float64; seed_py::Float64
    y::Float64; py::Float64
    prime::Int; T::Float64
    trace::Float64; detDT::Float64
    resnorm::Float64; iters::Int
    sec_y::Vector{Float64}; sec_py::Vector{Float64}
    history_y::Vector{Float64}; history_py::Vector{Float64}
    traj_x::Vector{Float64}; traj_y::Vector{Float64}
    class::String; index::Int
    sym_y::Vector{Float64}; sym_py::Vector{Float64}
    sym_stab::String; sym_mult::Int; sym_id::Int
    id::Int; origin::String
end

orbit_table() = DataFrame(Row[])

sec_matrix(o) = permutedims([o.sec_y o.sec_py])


# =====================================================================
#  Section geometry
# =====================================================================

"px^2 on the section x = 0, as a function of (y, py)."
px2(y, py, E, p) = 2p[2] * (E - HenonHeiles.potential(0.0, y, p)) - py^2

in_section(v, E, p; margin = 0.0) = px2(v[1], v[2], E, p) > margin

pymax(y, E, p) = sqrt(max(0.0, 2 * p[2] * (E - HenonHeiles.potential(0.0, y, p))))

"""
Lift (y, py) to a 4D state. The x-offset follows sign(px), so the trajectory
leaves the section immediately and no phantom crossing is recorded at t = 0.
Returns `nothing` outside the energy boundary.
"""
function lift(v, E, p; sgn = +1)
    a = px2(v[1], v[2], E, p)
    a <= 0 && return nothing
    return [sgn * EPS_OFF, v[1], sgn * sqrt(a), v[2]]
end


# =====================================================================
#  Integrators
# =====================================================================

"""
Build the fast (residuals) and dense (trajectories) integrators.

Returns the nmax Refs as well -- always construct `SectionParams` from these
returned Refs, never from freshly-made ones, or writes to `prm.nmax_*` will be
invisible to the callbacks.

Only upcrossings are recorded (affect_neg! = nothing), so every stored crossing
has px > 0. The symmetry bookkeeping relies on that.
"""
function create_integrators(p; tmax = 20_000.0, dt = DT, nfast = 1, ndense = 40)
    nmax_fast, nmax_dense = Ref(nfast), Ref(ndense)
    condition(u, t, integ) = u[1]

    yf, pyf, tsf = Float64[], Float64[], Float64[]
    affect_f!(integ) = begin
        push!(yf, integ.u[2]); push!(pyf, integ.u[4]); push!(tsf, integ.t)
        length(yf) ≥ nmax_fast[] && terminate!(integ)
    end
    integ_fast = init(ODEProblem(HenonHeiles.equations!, zeros(4), (0.0, tmax), p),
                      Vern9(); abstol = TOL_INTEGRATE, reltol = TOL_INTEGRATE,
                      save_everystep = false, save_start = false,
                      callback = ContinuousCallback(condition, affect_f!, nothing;
                                                    abstol = TOL_CALLBACK))

    yd, pyd, tsd = Float64[], Float64[], Float64[]
    affect_d!(integ) = begin
        push!(yd, integ.u[2]); push!(pyd, integ.u[4]); push!(tsd, integ.t)
        length(yd) ≥ nmax_dense[] && terminate!(integ)
    end
    integ_dense = init(ODEProblem(HenonHeiles.equations!, zeros(4), (0.0, tmax), p),
                       Vern9(); abstol = TOL_INTEGRATE, reltol = TOL_INTEGRATE, saveat = dt,
                       callback = ContinuousCallback(condition, affect_d!, nothing;
                                                     abstol = TOL_CALLBACK))

    return (; integ_fast, yf, pyf, tsf, nmax_fast,
              integ_dense, yd, pyd, tsd, nmax_dense)
end

function SectionParams(E, p; tmax = 20_000.0, dt = DT, nfast = 1, ndense = 40)
    b = create_integrators(p; tmax, dt, nfast, ndense)
    return SectionParams(b.integ_fast, b.yf, b.pyf, b.tsf, b.nmax_fast,
                         b.integ_dense, b.yd, b.pyd, b.tsd, b.nmax_dense,
                         E, p, tmax)
end

"Run body with prm.nmax_dense[] temporarily set to n, then restore it."
function with_ndense(f, prm, n)
    old = prm.nmax_dense[]
    prm.nmax_dense[] = n
    try
        return f()
    finally
        prm.nmax_dense[] = old
    end
end

"""
Integrate from a 4D state to the first section crossing that is safely inside
the energy boundary, and return it as (y, py). Crossings earlier than `tmin` are
skipped: a symmetry image can start arbitrarily close to the section.
"""
function flow_to_section(u0, prm; search = 12, margin = 1e-8, tmin = TMIN_CROSS)
    with_ndense(prm, search) do
        empty!(prm.yd); empty!(prm.pyd); empty!(prm.tsd)
        reinit!(prm.integ_dense, u0)
        solve!(prm.integ_dense)
        for i in eachindex(prm.yd)
            prm.tsd[i] > tmin || continue
            w = [prm.yd[i], prm.pyd[i]]
            in_section(w, prm.E, prm.p; margin) && return w
        end
        return nothing
    end
end


# =====================================================================
#  Poincare map and residual
# =====================================================================

"n-th crossing of the section, starting from v. This is T^n(v)."
function poincare_map(v, n::Int, prm::SectionParams)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")

    n > prm.nmax_fast[] && (prm.nmax_fast[] = n)

    empty!(prm.yf); empty!(prm.pyf); empty!(prm.tsf)
    reinit!(prm.integ_fast, u0)
    solve!(prm.integ_fast)

    length(prm.yf) < n &&
        error("only $(length(prm.yf)) crossings in t < $(prm.tmax) (need $n)")
    return [prm.yf[n], prm.pyf[n]]
end

Fres(v, n, prm) = poincare_map(v, n, prm) - v

"Finite everywhere: TrustRegion probes outside the boundary and NaN poisons it."
function Fres_safe(v, n, prm)
    a = px2(v[1], v[2], prm.E, prm.p)
    a <= 0 && return fill(1.0 + 100 * sqrt(-a), 2)
    return Fres(v, n, prm)
end

"Residual norm, or NaN if v is not numerically inside the section."
function residual_norm(v, n, prm)
    v === nothing && return NaN
    in_section(v, prm.E, prm.p) || return NaN
    try
        return norm(Fres(v, n, prm))
    catch
        return NaN
    end
end

"Dense integration from v. Returns (sol, section points as 2xN, crossing times)."
function section_trj(v, prm)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")

    empty!(prm.yd); empty!(prm.pyd); empty!(prm.tsd)
    reinit!(prm.integ_dense, u0)
    solve!(prm.integ_dense)
    return prm.integ_dense.sol, permutedims([prm.yd prm.pyd]), copy(prm.tsd)
end

"""
Smallest k with |T^k v - v| < tol.

Returns (; traj, pMap, Nperiod, Tperiod); Nperiod is `nothing` if the orbit did
not close within `search` crossings.
"""
function minPeriodicity(v, prm; tol = TOL_PERIOD, search = 40)
    with_ndense(prm, search) do
        sol, trace, ts = section_trj(v, prm)

        for i in axes(trace, 2)
            if norm(v .- trace[:, i]) < tol
                k = searchsortedfirst(sol.t, ts[i])
                return (; traj = sol.u[1:k], pMap = trace[:, 1:i],
                          Nperiod = i, Tperiod = ts[i])
            end
        end

        @warn "no closure within $search crossings" v tol
        return (; traj = sol.u, pMap = trace, Nperiod = nothing, Tperiod = nothing)
    end
end


# =====================================================================
#  Jacobian and root finding
# =====================================================================

function jacobian!(J, v, n, prm, d)
    E, p = prm.E, prm.p
    for (j, e) in enumerate(([d, 0.0], [0.0, d]))
        vp, vm   = v .+ e, v .- e
        okp, okm = in_section(vp, E, p), in_section(vm, E, p)
        if okp && okm
            J[:, j] = (Fres(vp, n, prm) .- Fres(vm, n, prm)) ./ (2d)   # central
        elseif okp
            J[:, j] = (Fres(vp, n, prm) .- Fres(v,  n, prm)) ./ d      # forward
        elseif okm
            J[:, j] = (Fres(v,  n, prm) .- Fres(vm, n, prm)) ./ d      # backward
        else
            error("both probes outside boundary at $v -- v is on the edge")
        end
    end
    return J
end

jacobian(v, n, prm, d = FD_STEP) = jacobian!(zeros(2, 2), v, n, prm, d)

get_DT(v, n, prm; d = FD_STEP) = I + jacobian(v, n, prm, d)

"""
Damped Newton with a step cap and backtracking line search. Kept as an
independent check on `solve_orbit`; unlike NonlinearSolve it records the
iterate history.
"""
function find_orbit(v0, n, prm;
                    N_max = 100, d = FD_STEP, tol = TOL_PO_ROOT,
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
            return (v = v, DT = DT, converged = true,
                    history = hist[:, 1:i], resnorm = rn,
                    comment = "|r| = $rn  det(DT) = $(det(DT))")
        end

        jacobian!(J, v, n, prm, d)
        κ = cond(J)
        κ > 1e12 && @warn "ill-conditioned Jacobian" iterate = i cond = κ

        step = J \ r
        sn   = norm(step)
        sn > dmax && (step .*= dmax / sn)       # cap: keeps Newton local

        λ, accepted = 1.0, false                # MUST start false
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
NonlinearSolve variant, used by the sweep. AutoFiniteDiff is mandatory: the ODE
callback rejects Duals. Returns an empty `history` -- use `find_orbit` if the
iterate path is needed.
"""
function solve_orbit(v0, n, prm; tol = TOL_PO_ROOT, maxiters = 300)
    v = collect(float.(v0))
    in_section(v, prm.E, prm.p) ||
        return (v = v, DT = nothing, converged = false, resnorm = Inf,
                history = zeros(2, 0), comment = "seed outside boundary")

    prob = NonlinearProblem((w, q) -> Fres_safe(w, n, q), v, prm)
    sol  = solve(prob, TrustRegion(autodiff = AutoFiniteDiff(fdtype = Val(:central)));
                 abstol = tol, maxiters)

    ok = SciMLBase.successful_retcode(sol) && in_section(sol.u, prm.E, prm.p)
    ok || return (v = sol.u, DT = nothing, converged = false,
                  resnorm = norm(sol.resid), history = zeros(2, 0),
                  comment = "$(sol.retcode)")

    DT = get_DT(sol.u, n, prm)
    return (v = sol.u, DT = DT, converged = true, resnorm = norm(sol.resid),
            history = zeros(2, 0), comment = "$(sol.retcode)  det(DT) = $(det(DT))")
end

"""
Solve from a seed and package the result as a `Row`. Returns `nothing` if
Newton fails or the orbit does not close.

The symmetry fields are left empty; `annotate_symmetry!` fills them once the
cheap raw dedup has thinned the table.
"""
function analyse_seed(v0, n, prm; id = 0, origin = "seed",
                      tol = TOL_PO_ROOT, search = 40)
    res = solve_orbit(v0, n, prm; tol)
    res.converged || return nothing

    orb = minPeriodicity(res.v, prm; search)
    orb.Nperiod === nothing && return nothing

    DT    = get_DT(res.v, orb.Nperiod, prm)
    τ     = tr(DT)
    index = kind_index(τ)

    return (; E = prm.E, n, seed_y = v0[1], seed_py = v0[2],
              y = res.v[1], py = res.v[2],
              prime = orb.Nperiod, T = orb.Tperiod,
              trace = τ, detDT = det(DT),
              resnorm = res.resnorm, iters = size(res.history, 2),
              sec_y = orb.pMap[1, :], sec_py = orb.pMap[2, :],
              history_y = res.history[1, :], history_py = res.history[2, :],
              traj_x = [u[1] for u in orb.traj], traj_y = [u[2] for u in orb.traj],
              class = KIND_LABEL[index], index,
              sym_y = Float64[], sym_py = Float64[],
              sym_stab = "", sym_mult = 0, sym_id = 0,
              id, origin)
end


# =====================================================================
#  Symmetry images of an orbit
# =====================================================================

"Is w within tol of any section point of the orbit?"
on_orbit(w, sec_y, sec_py; tol = TOL_ROOT_COMP) =
    any(i -> hypot(sec_y[i] - w[1], sec_py[i] - w[2]) < tol, eachindex(sec_y))

"""
    symmetry_data(v, sec_y, sec_py, prm)

Map the orbit through v by all 12 group elements and bring each image back to
the section.

Costs 5 integrations, not 11: the group splits into 6 cosets of the
section-stabilising subgroup {E, σΘ}, one representative of each needs a flow
back to the section, and the identity coset is free because v is already a root.
Within a coset the partner follows algebraically from (y, py) -> (y, -py).

Returns (; roots, stab, mult), where `roots` holds one section point per group
element with duplicates removed, `stab` names the elements that map the orbit to
itself, and `mult` = 12 / |stab| is the number of distinct orbits in the class.
`mult` is 0 when some image could not be brought back to the section, which
makes the count unreliable.
"""
function symmetry_data(v, sec_y, sec_py, prm;
                       search = 12, margin = 1e-8, tol = TOL_ROOT_COMP)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && return (roots = Vector{Float64}[], stab = String[], mult = 0)

    roots, stab, complete = Vector{Float64}[], String[], true

    for g in COSET_REPS
        w0 = g.name == "E" ? collect(float.(v)) :
             flow_to_section(apply_sym(g, u0), prm; search, margin)
        if w0 === nothing
            complete = false
            continue
        end
        for h in SEC_STAB
            w  = sec_action(h, w0)
            gh = compose(g, h)
            on_orbit(w, sec_y, sec_py; tol) && push!(stab, gh.name)
            any(r -> norm(r - w) < tol, roots) || push!(roots, w)
        end
    end

    mult = 0
    if complete && !isempty(stab)
        q, r = divrem(length(SYM_GROUP), length(stab))
        if r == 0
            mult = q
        else
            @warn "stabiliser order does not divide 12 -- tighten tol or margin" v stab
        end
    end
    return (; roots, stab, mult)
end

"""
Fill `sym_y`, `sym_py`, `sym_stab` and `sym_mult` for every row of `df`.

`df` must hold a single energy, matching `prm`. Run this only on a table that
`dedup_raw` has already thinned: it costs 5 dense integrations per row.
"""
function annotate_symmetry!(df, prm; tol = TOL_ROOT_COMP, search = 12, margin = 1e-8)
    for i in 1:nrow(df)
        s = try
            symmetry_data([df.y[i], df.py[i]], df.sec_y[i], df.sec_py[i], prm;
                          search, margin, tol)
        catch e
            @warn "symmetry annotation failed" id = df.id[i] exception = e
            continue
        end
        df.sym_y[i]    = [w[1] for w in s.roots]
        df.sym_py[i]   = [w[2] for w in s.roots]
        df.sym_stab[i] = join(s.stab, ",")
        df.sym_mult[i] = s.mult
    end
    return df
end

"Human-readable name for a stabiliser subgroup."
function sym_type(stab::AbstractString)
    els = Set(split(stab, ","; keepempty = false))
    isempty(els)      && return "unknown"
    length(els) == 12 && return "fully symmetric"
    has(x) = x in els
    rot = has("C3") && has("C3²")
    mir = has("σ")  || has("σC3")  || has("σC3²")
    rev = has("Θ")  || has("σΘ")   || has("C3Θ") || has("C3²Θ") ||
          has("σC3Θ") || has("σC3²Θ")
    rot && mir && return "C3v invariant"
    rot        && return "C3 invariant"
    mir && rev && return "mirror + reversible"
    mir        && return "mirror invariant"
    rev        && return "reversible"
    return "generic"
end


# =====================================================================
#  Deduplication
# =====================================================================

"""
Do two rows describe the same orbit up to a symmetry of the system?

The geometric test is decisive: does any symmetry image of either row land on
the other row's section points? `T` and `tr(DT)` are useful diagnostics, but
they are computed independently and must not reject a geometric match.
"""
function same_class(o1, o2; tol = TOL_ROOT_COMP)
    o1.prime == o2.prime                 || return false
    isapprox(o1.E, o2.E; atol = 1e-12)   || return false
    return symmetry_match(o1, o2; tol)
end

"Report invariant mismatches without using them to reject a geometric match."
function class_diagnostics(o1, o2; rtol_T = 1e-6, rtol_tr = 1e-4, atol_tr = 1e-6)
    return (; same_prime = o1.prime == o2.prime,
            same_energy = isapprox(o1.E, o2.E; atol = 1e-12),
            same_geometry = symmetry_match(o1, o2),
            same_period = isapprox(o1.T, o2.T; rtol = rtol_T),
            same_trace = isapprox(o1.trace, o2.trace;
                                  rtol = rtol_tr, atol = atol_tr))
end

"Check both directions because either row may have incomplete symmetry data."
function symmetry_match(o1, o2; tol = TOL_ROOT_COMP)
    same_point_set(o1, o2; tol) && return true
    any(a -> on_orbit([o1.sym_y[a], o1.sym_py[a]], o2.sec_y, o2.sec_py; tol),
        eachindex(o1.sym_y)) && return true
    return any(a -> on_orbit([o2.sym_y[a], o2.sym_py[a]], o1.sec_y, o1.sec_py; tol),
               eachindex(o2.sym_y))
end

"Literally the same orbit -- section points in common, no symmetry applied."
same_point_set(o1, o2; tol = TOL_ROOT_COMP) =
    any(a -> on_orbit([o1.sec_y[a], o1.sec_py[a]], o2.sec_y, o2.sec_py; tol),
        eachindex(o1.sec_y))

"""
Deduplicate under an arbitrary equivalence predicate by connected components.
Every pair is tested, so approximate equivalence cannot depend on row order.
`assign` receives the component representative and merged row indices.
"""
function dedup_by(df, equal; assign = nothing)
    n = nrow(df)
    parent = collect(1:n)

    function root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end

    for i in 2:n, j in 1:i-1
        equal(df[j, :], df[i, :]) || continue
        rj, ri = root(j), root(i)
        rj == ri && continue
        parent[ri] = rj
        assign === nothing || assign(rj, ri)
    end

    keep = falses(n)
    representatives = Dict{Int, Int}()
    for i in 1:n
        r = root(i)
        representatives[r] = get(representatives, r, i)
    end
    for i in values(representatives)
        keep[i] = true
    end
    return df[keep, :]
end

"Collapse literal duplicates within one energy. Cheap: no symmetry images needed."
dedup_raw(df; tol = TOL_ROOT_COMP) =
    dedup_by(df, (a, b) -> a.prime == b.prime && same_point_set(a, b; tol))

"""
Collapse a table to one row per symmetry class within one energy and stamp
`sym_id`. Rows are sorted by `prime` then `resnorm` first, so the shortest and
best-converged orbit becomes the representative.

Falls back to `dedup_raw` behaviour on rows that `annotate_symmetry!` has not
reached.
"""
function dedup_classes(df; tol = TOL_ROOT_COMP)
    isempty(df) && return copy(df)
    sorted = sort(DataFrame(df), [:prime, :resnorm])
    reps   = dedup_by(sorted, (a, b) -> same_class(a, b; tol))
    reps.sym_id .= 1:nrow(reps)
    return reps
end

"""
Full reduction of a multi-energy catalogue. Classes are only ever compared
within one energy, since section points at different E can coincide without
being the same orbit.
"""
function dedup_all(df; tol = TOL_ROOT_COMP, by_class = true)
    isempty(df) && return copy(df)
    f = by_class ? (d -> dedup_classes(d; tol)) :
                   (d -> dedup_raw(sort(DataFrame(d), [:prime, :resnorm]); tol))
    sort!(vcat([f(DataFrame(g)) for g in groupby(df, :E)]...), [:E, :prime, :T])
end

"Fast regression tests for deduplication; no ODE integrations are required."
function test_dedup_logic()
    row(; y, py, sec_y = [y], sec_py = [py], sym_y = Float64[],
        sym_py = Float64[], prime = 1, T = 2.0, trace = 0.5,
        resnorm = 1e-10) = (; 
        E = 0.1, n = prime, seed_y = y, seed_py = py,
        y = y, py = py, prime = prime, T = T, trace = trace,
        detDT = 1.0, resnorm = resnorm, iters = 1,
        sec_y = sec_y, sec_py = sec_py, history_y = [y], history_py = [py],
        traj_x = [0.0], traj_y = [y], class = "elliptic", index = 1,
        sym_y = sym_y, sym_py = sym_py, sym_stab = "", sym_mult = 0, sym_id = 0,
        id = 0, origin = "test")

    @testset "deduplication" begin
        @test length(SYM_GROUP) == 12
        @test length(SEC_STAB) == 2
        @test length(COSET_REPS) == 6

        literal_a = row(y = 0.1, py = 0.2)
        literal_b = row(y = 0.1 + 1e-8, py = 0.2 - 1e-8)
        @test same_point_set(literal_a, literal_b; tol = 1e-6)
        @test nrow(dedup_raw(DataFrame([literal_a, literal_b]); tol = 1e-6)) == 1

        partner_a = row(y = 0.3, py = 0.4,
                        sym_y = [0.7], sym_py = [-0.2])
        partner_b = row(y = 0.7, py = -0.2, T = 2.01, trace = 0.7)
        @test same_class(partner_a, partner_b; tol = 1e-6)
        @test same_class(partner_b, partner_a; tol = 1e-6)
        @test !class_diagnostics(partner_a, partner_b).same_period
        @test !class_diagnostics(partner_a, partner_b).same_trace

        # A symmetry annotation on only one row must still be sufficient.
        @test same_class(row(y = 0.7, py = -0.2), partner_a; tol = 1e-6)

        different_root = row(y = 0.8, py = -0.1)
        different_prime = row(y = 0.7, py = -0.2, prime = 2)
        different_energy = merge(partner_b, (; E = 0.2))
        @test !same_class(partner_a, different_root; tol = 1e-6)
        @test !same_class(partner_a, different_prime; tol = 1e-6)
        @test !same_class(partner_a, different_energy; tol = 1e-6)

        class_rows = DataFrame([
            partner_a,
            partner_b,
            row(y = -0.4, py = 0.2, resnorm = 1e-8),
        ])
        classes = dedup_classes(class_rows; tol = 1e-6)
        @test nrow(classes) == 2
        @test classes.sym_id == [1, 2]

        # Approximate matching can form a chain: A matches B and B matches C.
        chain = DataFrame([
            row(y = 0.0, py = 0.0),
            row(y = 0.1, py = 0.0),
            row(y = 0.2, py = 0.0),
        ])
        chain_dedup = dedup_by(chain, (a, b) -> abs(a.y - b.y) <= 0.11)
        @test nrow(chain_dedup) == 1
    end

    return true
end


# =====================================================================
#  Expanding a class back to its members (for plotting)
# =====================================================================

"""
The inverse of `dedup_classes`: regenerate every distinct member of each
symmetry class from its representative. Useful for drawing the full family; not
needed for analysis, where one orbit per class is the point.

Each stored image root is polished by a Newton solve -- it is already accurate,
so this converges in one or two steps -- and the results are deduplicated. A
mismatch against `sym_mult` is warned about rather than silently accepted.
"""
function expand_symmetry(df, prm; tol = TOL_ROOT_COMP, verbose = false)
    out = orbit_table()
    nid = 0
    for o in eachrow(df)
        members = orbit_table()
        for a in eachindex(o.sym_y)
            new = try
                analyse_seed([o.sym_y[a], o.sym_py[a]], o.prime, prm;
                             origin = "sym(id$(o.id))")
            catch e
                verbose && @warn "expansion failed" id = o.id exception = e
                nothing
            end
            new === nothing && continue
            new.prime == o.prime || continue
            push!(members, new)
        end
        members = dedup_raw(members; tol)
        members.sym_id   .= o.sym_id
        members.sym_mult .= o.sym_mult
        members.sym_stab .= o.sym_stab
        for k in 1:nrow(members)
            nid += 1
            members.id[k] = nid
        end
        o.sym_mult == 0 || nrow(members) == o.sym_mult ||
            @warn "expanded $(nrow(members)) members, expected $(o.sym_mult)" id = o.id
        append!(out, members)
    end
    return out
end


# =====================================================================
#  Energy grids
# =====================================================================

"""
    log_energies(Emin, Emax; n, per_decade, digits)

Geometric energy grid: a constant ratio ρ = E_{i+1}/E_i, so

    Δ(log E) = log ρ = const,

which is what makes the points evenly spaced on a log-E axis. Note that the
*absolute* step is not constant -- it scales as ΔE_i = (ρ-1) E_i, fine near
E = 0 and coarse near the top of the range.

Give either the total count `n` or the density `per_decade`; the latter is
usually the natural handle, since
`n = per_decade * log10(Emax/Emin) + 1`.

Returned in ascending order. `follow_orbits` sorts descending itself.
"""
function log_energies(Emin, Emax; n = nothing, per_decade = nothing)
    0 < Emin < Emax || throw(ArgumentError("need 0 < Emin < Emax"))
    (n === nothing) == (per_decade === nothing) &&
        throw(ArgumentError("give exactly one of n or per_decade"))
    decades = log10(Emax / Emin)
    n = n === nothing ? round(Int, per_decade * decades) + 1 : n
    if n < 2
        throw(ArgumentError(
            "per_decade = $per_decade over $(round(decades; sigdigits = 3)) decades " *
            "gives n = $n. The range [$Emin, $Emax] is too narrow for geometric " *
            "spacing in E to mean anything -- pass n explicitly, or use " *
            "crit_energies(Ec, dmin, dmax) if you are resolving the approach to a " *
            "critical energy."))
    end
    return collect(logrange(Emin, Emax, n))
end

"""
    crit_energies(Ec, dmin, dmax; kwargs...)

Geometric in the *distance to a critical energy*: Δ(log(Ec - E)) = const, so
the grid clusters towards Ec rather than towards zero. Use this when the
feature of interest is the divergence of periods as E -> Ec (for the standard
parameters a = m = ω = 1 the escape energy is Ec = 1/6).

Returned in ascending E.
"""
crit_energies(Ec, dmin, dmax; kwargs...) =
    sort(Ec .- log_energies(dmin, dmax; kwargs...))

"Report the sampling of an energy grid: ratio, log-step, and absolute extremes."
function grid_report(Es)
    E = sort(Es)
    ρ = E[2:end] ./ E[1:end-1]
    d = diff(log10.(E))
    @printf("%d energies, E in [%.4g, %.4g], %.2f decades\n",
            length(E), first(E), last(E), log10(last(E) / first(E)))
    @printf("  ratio  E_{i+1}/E_i : %.6f … %.6f\n", minimum(ρ), maximum(ρ))
    @printf("  Δlog10 E          : %.6f … %.6f  (constant iff geometric)\n",
            minimum(d), maximum(d))
    @printf("  ΔE                : %.3g … %.3g\n",
            minimum(diff(E)), maximum(diff(E)))
    return (; ratio = extrema(ρ), dlog = extrema(d))
end


# =====================================================================
#  Sweeps
# =====================================================================

"Seed grid on the section, inset from the energy boundary by `margin`."
function section_grid(E, p; ny = 10, npy = 10, margin = 0.03)
    roots = real.(filter(r -> abs(imag(r)) < 1e-10,
                         HenonHeiles.limit_of_initial_y0(E, p)))
    sort!(roots)
    ymin, ymax = roots[1], roots[2]      # the two turning points bounding the well
    dy = margin * (ymax - ymin)
    seeds = Vector{Float64}[]
    for y in range(ymin + dy, ymax - dy, length = ny)
        pm = pymax(y, E, p)
        pm <= 0 && continue
        for s in range(margin, 1 - margin, length = npy)
            push!(seeds, [y, s * pm])
        end
    end
    return seeds
end

function sweep!(df, seeds, n, prm; tol = TOL_PO_ROOT, verbose = false)
    nid = isempty(df) ? 0 : maximum(df.id)
    @showprogress dt = 1 desc = "E=$(round(prm.E; digits = 4)) n=$n" for v0 in seeds
        orb = try
            analyse_seed(v0, n, prm; id = nid + 1, tol)
        catch e
            verbose && @warn "seed failed" seed = v0 exception = e
            continue
        end
        orb === nothing && continue
        orb.prime == n  || continue
        nid += 1
        push!(df, merge(orb, (; id = nid)))
    end
    return df
end

"Seeds tagged with the prime period they were found at, as `follow_orbits` carries them."
const TaggedSeeds = Vector{Tuple{Vector{Float64}, Int}}

"Group tagged seeds into prime -> points."
function seeds_by_prime(pairs)
    by = Dict{Int, Vector{Vector{Float64}}}()
    for (v, k) in pairs
        push!(get!(by, Int(k), Vector{Vector{Float64}}()), collect(float.(v)))
    end
    return by
end

"""
    seed_router(seeds, primes, E, p; ny, npy, match)

Return a function `n -> points` that decides which seeds are tried when looking
for period-n orbits. Accepted forms of `seeds`:

  * `nothing`            -- build a fresh `section_grid`, used for every n
  * a flat list of `(y, py)` points -- used for every n
  * a list of `(point, prime)` pairs -- routed by prime
  * a flat list plus a parallel `primes` vector -- routed by prime

`match` controls the routing of tagged seeds:

  * `:exact`    -- only seeds whose prime is n. Cheapest, and correct for pure
                   continuation: a period-n orbit stays period-n as E varies.
  * `:divisors` -- seeds whose prime divides n. Costs more, but a period-n orbit
                   born in a period-multiplying bifurcation emerges next to its
                   period-n/k ancestor, so this is what catches new branches.
  * `:all`      -- every seed for every n, i.e. the old behaviour.
"""
function seed_router(seeds, primes, E, p; ny = 5, npy = 5, match = :exact)
    match in (:exact, :divisors, :all) ||
        throw(ArgumentError("match must be :exact, :divisors or :all"))

    if seeds === nothing
        grid = section_grid(E, p; ny, npy)
        return _ -> grid
    end

    tagged =
        if primes !== nothing
            length(primes) == length(seeds) ||
                throw(ArgumentError("seeds and primes must have equal length"))
            collect(zip(seeds, primes))
        elseif !isempty(seeds) && first(seeds) isa Tuple
            seeds
        else
            nothing
        end

    if tagged === nothing || match === :all
        flat = tagged === nothing ? collect(seeds) : [v for (v, _) in tagged]
        return _ -> flat
    end

    by = seeds_by_prime(tagged)
    empty_pts = Vector{Vector{Float64}}()

    match === :exact && return n -> get(by, n, empty_pts)
    return n -> reduce(vcat, (v for (k, v) in by if n % k == 0); init = empty_pts)
end

"""
Sweep all periods in `ns` at one energy, then reduce.

Seeds are routed by prime period (see `seed_router`), so a seed that was a
period-2 orbit at the previous energy is not re-solved against `F_1`, `F_3` and
`F_4`. With `ns = 1:4` and `match = :exact` this is a fourfold reduction in
Newton solves; the discarded work could only ever have produced rows that
`sweep!` rejects with `orb.prime == n`.

The reduction is staged so the expensive step runs on the fewest rows:
  1. `dedup_raw`          -- literal duplicates, pure section-point comparison
  2. `annotate_symmetry!` -- 5 integrations per surviving row
  3. `dedup_classes`      -- one row per symmetry class, `sym_id` stamped

`reduce_to` selects what comes back: `:classes` (default) returns the class
representatives, `:orbits` returns every distinct orbit with its symmetry data
filled in, `:raw` skips the reduction entirely.

Returns (; df, timing).
"""
function sweep_energy(E, ns, p; tmax = 5000.0, ny = 5, npy = 5,
                      ndense = 40, seeds = nothing, primes = nothing,
                      match = :exact, reduce_to = :classes, tol = TOL_ROOT_COMP)
    prm    = SectionParams(E, p; tmax, nfast = 1, ndense)
    pick   = seed_router(seeds, primes, E, p; ny, npy, match)
    df     = orbit_table()
    timing = DataFrame(E = Float64[], n = Int[], nseeds = Int[],
                       found = Int[], secs = Float64[])

    for n in ns
        pts = pick(n)
        if isempty(pts)
            push!(timing, (E, n, 0, 0, 0.0))     # keep the timing table rectangular
            continue
        end
        before = nrow(df)
        secs   = @elapsed sweep!(df, pts, n, prm)
        push!(timing, (E, n, length(pts), nrow(df) - before, secs))
    end

    reduce_to === :raw && return (; df, timing)

    df = dedup_raw(df; tol)
    annotate_symmetry!(df, prm; tol)
    reduce_to === :classes && (df = dedup_classes(df; tol))

    return (; df, timing)
end

"""
Run `sweep_energy` across energies on all available threads.

One SectionParams per task, created inside the task -- never shared, because
the integrators and their callback buffers are mutable state.
"""
function run_sweep_mt(Es, ns, p; tmax = 5000.0, ny = 5, npy = 5,
                      ndense = 40, reduce_to = :classes)
    println("running sweep, parallel computation\nreduce_to = $reduce_to\n")

    tasks = map(Es) do E
        Threads.@spawn sweep_energy(E, ns, p; tmax, ny, npy, ndense, reduce_to)
    end
    results = fetch.(tasks)

    df     = reduce(vcat, r.df     for r in results)
    timing = reduce(vcat, r.timing for r in results)

    sort!(df, [:E, :prime, :T])
    return (; df, timing)
end

"Tagged seed list from a catalogue: one (point, prime) pair per row."
seeds_from(d)::TaggedSeeds = [([d.y[i], d.py[i]], d.prime[i]) for i in 1:nrow(d)]

"""
Continuation: sweep the highest energy densely, then carry the orbits found at
each energy down as seeds for the next.

Two reductions compound here. Only one representative per symmetry class is
followed, shrinking the seed list by up to a factor of 12; and each seed is
tried only against its own prime period, saving a further factor of |ns|. The
symmetric partners are recoverable at any energy with `expand_symmetry`.

`match = :exact` is the fast path and assumes the branch structure is already
captured by the initial sweep. Use `match = :divisors` if you expect new
branches to appear as E decreases -- a period-n orbit born in a bifurcation sits
next to its period-n/k ancestor, and `:exact` will not look there.

If an energy yields nothing the previous seeds are reused rather than dropped,
so a single failed step does not end the continuation.

Returns (; df, timing).
"""
function follow_orbits(Es, ns, p; tmax = 5000.0, ny_init = 15, npy_init = 10,
                       ndense = 40, match = :exact, reduce_to = :classes,
                       tol = TOL_ROOT_COMP, verbose = true)
    Es = sort(Es; rev = true)            # do not mutate the caller's vector
    println("Following periodic orbits from E=$(Es[1]) to $(Es[end])")
    println("seeding with symmetry-class representatives " *
            "(reduce_to = $reduce_to, match = $match)\n")

    first_run = sweep_energy(Es[1], ns, p; tmax, ny = ny_init, npy = npy_init,
                             ndense, reduce_to, tol)
    seeds = seeds_from(first_run.df)

    all_df     = [first_run.df]
    all_timing = [first_run.timing]

    for e in Es[2:end]
        r = sweep_energy(e, ns, p; tmax, ndense, seeds, match, reduce_to, tol)
        push!(all_df, r.df)
        push!(all_timing, r.timing)
        if nrow(r.df) == 0
            verbose && @warn "no orbits survived at E = $e; keeping previous seeds"
        else
            seeds = seeds_from(r.df)
        end
    end

    df     = reduce(vcat, all_df)
    timing = reduce(vcat, all_timing)
    sort!(df, [:E, :prime, :T])
    sort!(timing, [:E])
    return (; df, timing)
end


# =====================================================================
#  Reporting
# =====================================================================

function report(df)
    header() = begin
        @printf("%3s %4s %5s %9s %9s %10s %10s  %-18s %4s %4s  %-20s %s\n",
                "id", "n", "prm", "y", "py", "T", "tr(DT)",
                "class", "sym", "mult", "type", "origin")
        println("-"^124)
    end

    println("="^124); header()
    for o in eachrow(df)
        @printf("%3d %4d %5d %9.5f %9.5f %10.4f %10.4f  %-18s %4d %4d  %-20s %s\n",
                o.id, o.n, o.prime, o.y, o.py, o.T, o.trace,
                o.class, o.sym_id, o.sym_mult, sym_type(o.sym_stab), o.origin)
    end
    println("="^124); header()

    bad = filter(o -> abs(o.detDT - 1) > 1e-4, eachrow(df))
    isempty(bad) ||
        @warn "rows with det(DT) far from 1 -- check the FD step" ids = [o.id for o in bad]

    miss = filter(o -> o.sym_mult == 0, eachrow(df))
    isempty(miss) ||
        @warn "rows with an incomplete symmetry class" ids = [o.id for o in miss]

    if !isempty(df) && any(df.sym_mult .> 0)
        @printf("\n%d representatives  ->  %d orbits after expansion\n",
                nrow(df), sum(max(o.sym_mult, 1) for o in eachrow(df)))
    end
end

"Count of symmetry types present, per energy and prime period."
sym_summary(df) = combine(groupby(transform(df, :sym_stab => ByRow(sym_type) => :type),
                                  [:E, :prime, :type]), nrow => :count)


# =====================================================================
#  Selecting and linking rows
# =====================================================================

"""
    at_energy(df, E; rtol)

Rows of a multi-energy catalogue at one energy. Every single-energy plot goes
through this: the section is a different manifold at every E, so drawing rows
from several energies on one set of axes superimposes an orbit on itself once
per energy and looks like a deduplication failure.
"""
function at_energy(df, E; rtol = 1e-8)
    sub = filter(row -> isapprox(row.E, E; rtol), df)
    if nrow(sub) == 0 && nrow(df) > 0
        Es = unique(df.E)
        error("no rows at E = $E; the table holds $(length(Es)) energies " *
              "in [$(minimum(Es)), $(maximum(Es))]")
    end
    return sub
end

"""
    trace_branches(df; tol, rtol_T)

Link rows across energies into continuation branches and return a copy with a
`:branch` column.

`sym_id` cannot do this job: `dedup_classes` assigns it independently at each
energy as `1:nrow`, so the same physical branch can be `sym_id = 3` at one
energy and `7` at the next. Branches are instead traced geometrically, walking
energies in the continuation direction (high to low) and matching each row to
its nearest predecessor of the same prime period.

The match is symmetry-aware: the distance to a predecessor is minimised over
that predecessor's stored class roots, so a branch is still followed when
`dedup_classes` happens to pick a different member of the symmetry class as
representative at the next energy.
"""
function trace_branches(df; tol = 0.05, rtol_T = 0.2)
    d = sort(DataFrame(df), [:E, :prime, :T])
    d.branch = zeros(Int, nrow(d))
    isempty(d) && return d

    "distance from row i's root to any class root of row j"
    function dist(i, j)
        v = (d.y[i], d.py[i])
        if isempty(d.sym_y[j])
            return hypot(d.y[j] - v[1], d.py[j] - v[2])
        end
        minimum(a -> hypot(d.sym_y[j][a] - v[1], d.sym_py[j][a] - v[2]),
                eachindex(d.sym_y[j]))
    end

    nb   = 0
    prev = Int[]
    for E in reverse(unique(d.E))          # continuation ran downwards in E
        cur  = findall(==(E), d.E)
        used = Set{Int}()
        for i in cur
            best, bd = 0, Inf
            for j in prev
                (d.prime[j] == d.prime[i] && j ∉ used)      || continue
                isapprox(d.T[i], d.T[j]; rtol = rtol_T)     || continue
                dj = dist(i, j)
                dj < bd && ((best, bd) = (j, dj))
            end
            if best > 0 && bd < tol
                d.branch[i] = d.branch[best]
                push!(used, best)
            else
                nb += 1
                d.branch[i] = nb
            end
        end
        prev = cur
    end
    return d
end

"How many energies each branch survives -- short branches are usually mismatches."
branch_report(d) = sort(combine(groupby(d, [:branch, :prime]),
                                nrow => :nE,
                                :E => minimum => :Emin,
                                :E => maximum => :Emax),
                        [:prime, :nE], rev = [false, true])


# =====================================================================
#  Plotting
# =====================================================================

"Load the background section scatter from a long simulation, or `nothing`."
function load_background(file = DATA_FILE)
    isfile(file) || (@warn "no background file" file; return nothing)
    y, py = Float64[], Float64[]
    for d in load(file, "results")
        append!(y,  d.sec_y)
        append!(py, d.sec_py)
    end
    return (y, py)
end

function config_backdrop!(ax, E, p; r = RGRID, levels = LEVELS)
    epot = [HenonHeiles.potential(x, y, p) for x in r, y in r]
    contour!(ax, r, r, epot; levels, colormap = :hsv, labels = true, linewidth = 0.8)
    contour!(ax, r, r, epot; levels = [E], color = C_CREAM, linewidth = 2.5)
    return ax
end

function section_backdrop!(ax, E, p; bg = nothing, seeds = nothing, nb = 200)
    bg === nothing || scatter!(ax, bg[1], bg[2]; color = (:grey, 0.45), markersize = 1.5)
    yb, pb = HenonHeiles.section_boundary_ranges(E, p, nb)
    scatter!(ax, HenonHeiles.section_boundary(yb, pb); color = C_CREAM, markersize = 3)
    seeds === nothing || scatter!(ax, first.(seeds), last.(seeds);
                                  color = (:white, 0.25), markersize = 4)
    return ax
end

"Legend for the categorical stability channel."
function stability_legend!(pos, kinds; lines = true, markers = true)
    elems = map(kinds) do k
        e = Any[]
        lines   && push!(e, LineElement(linestyle = KIND_LS[k], linewidth = 2, color = C_CREAM))
        markers && push!(e, MarkerElement(marker = KIND_MS[k], markersize = 12, color = C_CREAM))
        length(e) == 1 ? only(e) : e
    end
    Legend(pos, elems, KIND_LABEL[kinds], "stability";
           orientation = :horizontal, framevisible = false, tellheight = true)
end

function plot_config(df, prm; title = "configuration space")
    df  = at_energy(df, prm.E)
    fig = Figure(size = (1400, 900))
    ax  = Axis(fig[1, 1], xlabel = "x", ylabel = "y",
               title = title, aspect = DataAspect())
    config_backdrop!(ax, prm.E, prm.p)

    for o in eachrow(df)
        lines!(ax, o.traj_x, o.traj_y; color = pick_color(o.id),
               linestyle = KIND_LS[o.index],
               label = "id$(o.id) n=$(o.prime) T=$(round(o.T; digits = 1))")
    end
    stability_legend!(fig[2, 1], sort(unique(df.index)); markers = false)
    return fig, ax
end

function plot_section(df, prm; seeds = nothing, bg = nothing,
                      title = "surface of section")
    df  = at_energy(df, prm.E)
    fig = Figure(size = (1400, 900))
    ax  = Axis(fig[1, 1], xlabel = "y", ylabel = "p_y", title = title)
    section_backdrop!(ax, prm.E, prm.p; bg, seeds)

    for o in eachrow(df)
        scatter!(ax, o.sec_y, o.sec_py; color = pick_color(o.id), markersize = 12,
                 marker = KIND_MS[o.index],
                 label = "id$(o.id) n=$(o.prime) $(o.class)")
    end
    stability_legend!(fig[2, 1], sort(unique(df.index)); lines = false)
    return fig, ax
end

"""
    plot_orbits(df, prm; seeds, bg, cmap, label_ids, click_tol)

Two linked views of the catalogue: configuration space and the surface of
section. Returns `(; fig, ax1, ax2, alphas, info)` so the selection can also be
driven from the REPL, e.g. `o.alphas[3][] = 0.08` or
`foreach(a -> a[] = 1.0, o.alphas)`.

Pass class representatives to see one orbit per type, or the output of
`expand_symmetry` to see the full families.
"""
function plot_orbits(df, prm;
                     seeds     = nothing,
                     bg        = nothing,
                     cmap      = :viridis,
                     label_ids = true,
                     click_tol = 0.02,
                     verbose   = true)

    isempty(df) && error("nothing to plot: the orbit table is empty")

    # A section plot is a single-energy object. Without this the whole
    # continuation is drawn on one pair of axes and every orbit appears once
    # per energy, a hair apart -- which reads as a deduplication failure.
    nE = length(unique(df.E))
    df = at_energy(df, prm.E)
    verbose && nE > 1 &&
        @info "plot_orbits: $(nrow(df)) rows at E = $(prm.E), of $nE energies in the table"

    fig = Figure(size = (1800, 1000))
    ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = "y",
               title = "configuration space", aspect = DataAspect())
    ax2 = Axis(fig[1, 2], xlabel = "y", ylabel = "p_y",
               title = "surface of section  (E = $(round(prm.E; digits = 4)))")

    # left-drag is rectangle zoom by default and would swallow our clicks
    deregister_interaction!(ax1, :rectanglezoom)
    deregister_interaction!(ax2, :rectanglezoom)

    config_backdrop!(ax1, prm.E, prm.p; r = range(-1.0, 1.0, length = 220))
    section_backdrop!(ax2, prm.E, prm.p; bg, seeds)

    lo, hi = extrema(df.T)
    hi ≈ lo && (hi = lo + 1)                       # single orbit / degenerate range
    tcol(T) = get(colorschemes[cmap], (T - lo) / (hi - lo))

    alphas = [Observable(1.0) for _ in 1:nrow(df)]

    for (j, o) in enumerate(eachrow(df))
        b   = tcol(o.T)
        col = @lift(RGBAf(b.r, b.g, b.b, $(alphas[j])))

        lines!(ax1, o.traj_x, o.traj_y;
               color = col, linestyle = KIND_LS[o.index], linewidth = 2)

        scatter!(ax2, o.sec_y, o.sec_py;
                 color = col, marker = KIND_MS[o.index], markersize = 14,
                 strokewidth = 0.5, strokecolor = (:black, 0.6))

        if label_ids
            i0 = argmax(o.sec_py)                  # topmost point: labels spread out
            text!(ax2, o.sec_y[i0], o.sec_py[i0];
                  text = string(o.id), color = col, fontsize = 13,
                  offset = (8, 8), align = (:left, :bottom))
        end
    end

    Colorbar(fig[1, 3]; limits = (lo, hi), colormap = cmap, label = "orbit period T")
    stability_legend!(fig[2, 1:2], sort(unique(df.index)))

    info = Observable("click a trajectory or a section point to select an orbit")
    Label(fig[3, 1:3], info; tellwidth = false, fontsize = 15, font = :regular)

    btn = Button(fig[2, 3], label = "show all")
    on(_ -> foreach(a -> a[] = 1.0, alphas), btn.clicks)

    toggle!(j) = (alphas[j][] = alphas[j][] > 0.5 ? 0.08 : 1.0)

    function nearest(pos, xs, ys)
        best, bestd = 0, Inf
        for (j, o) in enumerate(eachrow(df))
            X, Y = xs(o), ys(o)
            for i in eachindex(X)
                d = hypot(X[i] - pos[1], Y[i] - pos[2])
                d < bestd && ((best, bestd) = (j, d))
            end
        end
        return best, bestd
    end

    function attach(ax, xs, ys)
        on(events(ax.scene).mousebutton) do ev
            (ev.button == Mouse.left && ev.action == Mouse.press) || return
            pos  = mouseposition(ax.scene)
            tol  = click_tol * maximum(widths(ax.finallimits[]))   # scales with zoom
            j, d = nearest(pos, xs, ys)
            (j > 0 && d < tol) || return
            toggle!(j)
            o = df[j, :]
            info[] = "id $(o.id)   prime = $(o.prime) (n = $(o.n))   " *
                     "T = $(round(o.T; digits = 4))   tr = $(round(o.trace; digits = 4))   " *
                     "$(o.class)   sym $(o.sym_id) x$(o.sym_mult) $(sym_type(o.sym_stab))   " *
                     "v = ($(round(o.y; digits = 6)), $(round(o.py; digits = 6)))   " *
                     "|r| = $(o.resnorm)   det = $(round(o.detDT; digits = 6))   [$(o.origin)]"
        end
    end

    attach(ax1, o -> o.traj_x, o -> o.traj_y)      # click a trajectory
    attach(ax2, o -> o.sec_y,  o -> o.sec_py)      # click a section point

    rowsize!(fig.layout, 1, Relative(0.85))
    return (; fig, ax1, ax2, alphas, info)
end

"All orbits of one prime period at one energy, in both views."
function plot_period(orbits, k; E = 0.1127, p = PARAM, bg = nothing, cmap = :managua,
                     show_newton = false, save_fig = false,
                     save_dir = FIG_DIR, show = true)
    sub = filter(row -> row.E ≈ E && row.prime == k, orbits)
    nrow(sub) == 0 && return nothing
    cols = resample_cmap(cmap, max(nrow(sub), 2))

    fc  = Figure(size = (1400, 900))
    axc = Axis(fc[1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect(),
               title = "config space — prime = $k, $(nrow(sub)) orbits, E = $E")
    config_backdrop!(axc, E, p)
    for (i, row) in enumerate(eachrow(sub))
        lines!(axc, row.traj_x, row.traj_y;
               color = cols[i], linewidth = 1.2, linestyle = KIND_LS[row.index])
    end
    stability_legend!(fc[2, 1], sort(unique(sub.index)); markers = false)

    lo, hi = extrema(sub.T)
    hi ≈ lo && (hi = lo + max(1e-9, abs(lo) * 1e-6))
    Colorbar(fc[1, 2]; limits = (lo, hi), colormap = cmap, label = "T")

    fp  = Figure(size = (1400, 900))
    axp = Axis(fp[1, 1], xlabel = "y", ylabel = "pᵧ", title = "section — prime = $k")
    section_backdrop!(axp, E, p; bg)
    for (i, row) in enumerate(eachrow(sub))
        show_newton && !isempty(row.history_y) &&
            lines!(axp, row.history_y, row.history_py;
                   color = (C_GOLD, 0.4), linewidth = 0.8)
        scatter!(axp, row.sec_y, row.sec_py;
                 color = cols[i], markersize = 9, strokewidth = 0.5,
                 marker = KIND_MS[row.index])
    end
    stability_legend!(fp[2, 1], sort(unique(sub.index)); lines = false)

    if save_fig
        mkpath(save_dir)
        tag = @sprintf("E%06.4f_p%02d", E, k)
        save(joinpath(save_dir, "config_$tag.png"),  fc; px_per_unit = 2)
        save(joinpath(save_dir, "section_$tag.png"), fp; px_per_unit = 2)
    end
    show && (display(GLMakie.Screen(), fc); display(GLMakie.Screen(), fp))
    return (; fc, axc, fp, axp, n = nrow(sub))
end

"All prime periods at one energy on a single section plot."
function plot_section_all(orbits, E; p = PARAM, bg = nothing, ks = nothing,
                          cmap = :managua, label = true,
                          save_fig = false, save_dir = FIG_DIR, show = true)
    sub = filter(row -> row.E ≈ E, orbits)
    nrow(sub) == 0 && return nothing
    ks   = something(ks, sort(unique(sub.prime)))
    sub  = filter(row -> row.prime in ks, sub)
    cols = resample_cmap(cmap, max(length(ks), 2))
    cidx = Dict(k => i for (i, k) in enumerate(ks))

    fp  = Figure(size = (1400, 900))
    axp = Axis(fp[1, 1], xlabel = "y", ylabel = "pᵧ",
               title = "section — E = $(round(E, digits = 4)), $(nrow(sub)) orbits")
    section_backdrop!(axp, E, p; bg)

    for row in eachrow(sub)
        c = cols[cidx[row.prime]]
        scatter!(axp, row.sec_y, row.sec_py; color = c, marker = KIND_MS[row.index],
                 markersize = 11, strokewidth = 0.5, strokecolor = :black)
        label && text!(axp, row.sec_y, row.sec_py;
                       text = fill(string(row.prime), length(row.sec_y)),
                       color = c, fontsize = 11, align = (:left, :bottom),
                       offset = (6, 4))
    end

    stability_legend!(fp[2, 1], sort(unique(sub.index)); lines = false)
    Colorbar(fp[1, 2]; limits = (minimum(ks) - 0.5, maximum(ks) + 0.5),
             colormap = cgrad(cmap, length(ks); categorical = true),
             ticks = ks, label = "prime period")

    if save_fig
        mkpath(save_dir)
        save(joinpath(save_dir, @sprintf("sectionall_E%06.4f.png", E)), fp; px_per_unit = 2)
    end
    show && display(GLMakie.Screen(), fp)
    return (; fp, axp, n = nrow(sub))
end

"""
    section_slider(orbits, p; ...)

Scan the surface of section across energies with a slider. One scatter series
per stability class, so the number of orbits can vary between energies.
"""
function section_slider(orbits, p; cmap = :viridis, markersize = 12)
    Es = sort(unique(orbits.E))
    isempty(Es) && error("no energies in the table")

    # ---- precompute: for each energy index, one (points, colours) pair per class
    npts  = [[Point2f[] for _ in 1:4] for _ in eachindex(Es)]
    ncols = [[Float64[] for _ in 1:4] for _ in eachindex(Es)]
    ninfo = Vector{String}(undef, length(Es))
    bnds  = Vector{Vector{Point2f}}(undef, length(Es))
    nlbl  = [String[]  for _ in eachindex(Es)]
    nlpos = [Point2f[] for _ in eachindex(Es)]

    Eidx = Dict(E => i for (i, E) in enumerate(Es))

    for r in eachrow(orbits)
        i, k = Eidx[r.E], r.index
        for (y, py) in zip(r.sec_y, r.sec_py)
            push!(npts[i][k],  Point2f(y, py))
            push!(ncols[i][k], r.T)
        end
        j = argmax(r.sec_py)                       # topmost point spreads labels out
        push!(nlpos[i], Point2f(r.sec_y[j], r.sec_py[j]))
        push!(nlbl[i],  string(r.prime))
    end

    for (i, E) in enumerate(Es)
        yb, pb  = HenonHeiles.section_boundary_ranges(E, p, 200)
        bnds[i] = HenonHeiles.section_boundary(yb, pb)
        s = filter(:E => ==(E), orbits)
        ninfo[i] = "E = $(round(E; digits = 5))   —   $(nrow(s)) orbits, " *
                   "periods $(isempty(s) ? "-" : "$(minimum(s.prime))–$(maximum(s.prime))")"
    end

    # ---- plot
    fig = Figure(size = (1500, 1000))
    sl  = Slider(fig[2, 1], range = eachindex(Es), startvalue = 1)
    idx = sl.value
    title_obs = @lift(ninfo[$idx])
    ax  = Axis(fig[1, 1], xlabel = "y", ylabel = "p_y", title = title_obs)


    scatter!(ax, @lift(bnds[$idx]); color = C_CREAM, markersize = 3)

    crange = extrema(orbits.T)
    for k in 1:4
        scatter!(ax, @lift(npts[$idx][k]);
                 color = @lift(ncols[$idx][k]),
                 colormap = cmap, colorrange = crange,
                 marker = KIND_MS[k], markersize,
                 strokewidth = 0.5, strokecolor = (:black, 0.5),
                 label = KIND_LABEL[k])
        text!(ax, @lift(nlpos[$idx]); text = @lift(nlbl[$idx]),
                color = C_CREAM, fontsize = 12,
                align = (:left, :bottom), offset = (6, 4))
    end

    
    # ax = Axis(fig[1, 1], xlabel = "y", ylabel = "p_y", title = title_obs)
    Colorbar(fig[1, 2]; limits = crange, colormap = cmap, label = "T")

    ymax = maximum(abs, Iterators.flatten(orbits.sec_y))
    pmax = maximum(abs, Iterators.flatten(orbits.sec_py))
    limits!(ax, -1.1ymax, 1.1ymax, -1.1pmax, 1.1pmax)
    limits!(ax, -1.1ymax, 1.1ymax, -1.1pmax, 1.1pmax)
    Legend(fig[1, 3],
           [MarkerElement(marker = KIND_MS[k], markersize = 12, color = :white)
            for k in 1:4],
           KIND_LABEL, "stability";
           framevisible = false)
    on(events(fig.scene).keyboardbutton) do ev
        ev.action == Keyboard.press || return
        i = sl.value[]
        ev.key == Keyboard.right && set_close_to!(sl, min(i + 1, length(Es)))
        ev.key == Keyboard.left  && set_close_to!(sl, max(i - 1, 1))
    end

    return (; fig, ax, slider = sl, Es)
end


"""
    plot_T_vs_E(orbits; ks, xlog, ylog, ...)

Period against energy, one line per continuation branch. This is what the
geometric energy grid from `log_energies` is for: with Δ(log E) constant the
points are evenly spaced along the x-axis at every scale.

Branches come from `trace_branches`, not from `sym_id` -- the latter is
reassigned independently at each energy and would scramble the curves. Pass a
table that already has a `:branch` column to skip the tracing.

Marker shape carries the stability class, so a branch switching from `:circle`
to `:xcross` marks the energy at which it goes unstable.
"""
function plot_T_vs_E(orbits; ks = nothing, xlog = false, ylog = true,
                     cmap = :managua, markers = true, min_len = 1,
                     tol = 0.05, save_fig = true, save_dir = FIG_DIR, show = true)
    isempty(orbits) && return nothing
    df = ks === nothing ? orbits : filter(row -> row.prime in ks, orbits)
    nrow(df) == 0 && return nothing

    "branch" in names(df) || (df = trace_branches(df; tol))

    branches = [b for b in groupby(sort(df, [:branch, :E]), :branch)
                if nrow(b) >= min_len]
    isempty(branches) && return nothing

    primes = sort(unique(df.prime))
    cols   = resample_cmap(cmap, max(length(primes), 2))
    cidx   = Dict(k => i for (i, k) in enumerate(primes))

    fig = Figure(size = (1400, 900))
    ax  = Axis(fig[1, 1];
               xlabel = "T", ylabel = "E",
               xscale = xlog ? log10 : identity,
               yscale = ylog ? log10 : identity,
               title  = "prime period vs energy ($(length(branches)) branches, " *
                        "$(length(unique(df.E))) energies)")

    for b in branches
        c = cols[cidx[first(b.prime)]]
        nrow(b) > 1 && lines!(ax, b.T, b.E,; color = c, linewidth = 1.4)
        markers && for r in eachrow(b)
            scatter!(ax, [r.E], [r.T]; color = c, markersize = 6,
                     marker = KIND_MS[r.index])
        end
    end

    stability_legend!(fig[2, 1], sort(unique(df.index)); lines = false)
    Colorbar(fig[1, 2]; limits = (minimum(primes) - 0.5, maximum(primes) + 0.5),
             colormap = cgrad(cmap, length(primes); categorical = true),
             ticks = primes, label = "prime period n")

    if save_fig
        mkpath(save_dir)
        save(joinpath(save_dir, "T_vs_E.png"), fig; px_per_unit = 2)
    end
    show && display(GLMakie.Screen(), fig)
    return (; fig, ax, df, nbranches = length(branches))
end


# =====================================================================
#  Driver
# =====================================================================

const NMAX_SEARCH = 20          # crossings the dense integrator may take

function main(; Emax = 0.166666, Emin = 1e-2, n_energies = 166*10, ns = 1:4,
                tmax = 10_000.0, match = :exact, reduce_to = :classes,
                outfile = "following_orbits_classes.jld2", save_data=true)

    Es = collect(range(Emin, Emax; length = n_energies))  
    println()

    res = follow_orbits(Es, ns, PARAM; tmax, ny_init = 15, npy_init = 10,
                        ndense = NMAX_SEARCH, match, reduce_to)

    orbits, timing = res.df, res.timing
    report(orbits)

    mkpath(SAVE_DATA_DIR)
    println("saving data: $save_data, file: $(joinpath(SAVE_DATA_DIR, outfile))")
    save_data && jldsave(joinpath(SAVE_DATA_DIR, outfile);
            orbits, timing, energies = Es,
            reduce_to = String(reduce_to), match = String(match),
            group_order = length(SYM_GROUP),
            newton_tol = TOL_PO_ROOT, ode_abstol = TOL_INTEGRATE, param = PARAM)

    return (; orbits, timing, Es)
end

# test_dedup_logic()

# res = main(
#     Emax = 0.166666,
#     Emin = 1e-2,
#     n_energies = 166*100,
#     ns = 1:9,
#     tmax = 20_000.0,
#     match = :exact,
#     reduce_to = :classes,
#     outfile = "05Sept_test_following_orbits_classes.jld2",
#     save_data = true
# )

# prm = SectionParams(Emax, PARAM;
#                     tmax = 20_000.0,
#                     ndense = NMAX_SEARCH)

# res_high = at_energy(res.orbits, prm.E)

# full = expand_symmetry(res_high, prm;
#                        tol = TOL_ROOT_COMP,
#                        verbose = true)

# plt_full = plot_orbits(full, prm; bg = nothing)
# display(plt_full.fig)

# plt_res = plot_orbits(res_high, prm; bg = nothing)
# display(plt_res.fig)