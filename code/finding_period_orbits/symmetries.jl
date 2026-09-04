import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2,
      NonlinearSolve, ADTypes, DataFrames, Dates, ProgressMeter, Printf
include("../styles/makie_theme.jl")
set_theme!(QUARTO_THEME)


const CONFIG_DIR = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
const DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
const FIG_DIR    = joinpath(@__DIR__, "../../figures/henon-heiles/periodic-orbits/")
data_file        = joinpath(DATA_DIR, "E0.1127-T10000.0-py0.0-n256.jld2")


# tolerances, precissions for 
const EPS_OFF       = 1e-10   # offset of u[1] to the surface of section, less than tol on Root
const TOL_INTEGRATE = 1e-14   # integrator precission
const TOL_CALLBACK  = 1e-13   # callback precission, less than integrator precission (!)
const TOL_PO_ROOT   = 1e-11   # closeness of to find the periodic orbits, i.e. T(v)-v < tol
const TOL_ROOT_COMP = 1e-7    # when comparing roots of PMAP to see if they are 'the same' vi-vj< tol


"returns -> 1 := elliptic, 2 := hyperbollic (stable), 3 := hyperbollic (unstable), 4 := parabolic,"
function kind_index(τ; ε = 1e-6)
    abs(abs(τ) - 2) < ε && return 4   # parabolic
    abs(τ) < 2 && return 1            # elliptic
    τ > 0 ? 2 : 3                     # hyperbollic (unstable, stable)
end

const KIND_LS    = [:solid, :dash, :dashdot, :dot]
const KIND_MS    = [:circle, :xcross, :diamond, :utriangle]
const KIND_LABEL = ["elliptic", "hyperbolic", "inverse hyperbolic", "parabolic"]






# Data Structures

struct SimParams{IF, ID, P}
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



const Row = @NamedTuple begin
    E::Float64; n::Int
    prime::Int; T::Float64
    trace::Float64; detDT::Float64
    class::String; index::Int
    resnorm::Float64; 
    sec_traj::Matrix{Float64}
    state_traj::Vector{Vector{Float64}}
    id::Int; origin::String
end
 
""" 
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
    id::Int; origin::String
"""

orbit_table() = DataFrame(Row[])


Tkin(p_i; p=(1.0, 1.0,1.0)) = p_i^2/(2*p[2])

Vpot(x,y; p=(1.0,1.0,1.0))  = p[2]*p[3]^2/2 *(x^2+y^2) + p[1]*(x^2*y - y^3/3)


function px2(v, E, p)
    x   = EPS_OFF
    y   = v[1]
    py  = v[2]
    m   = p[2]
    return (E-Tkin(py; p)-Vpot(x,y ;p)) * (2*m)
end 

function pymax(v, E, p)
    y   = v[1]
    m   = p[2]
    return sqrt( max(0.0, E - Vpot(EPS_OFF,y; p))) * 2*m
end

function on_section(v,E,p)
    py_test = v[2]
    py_on_sec= (py_test > pymax(v,E,p))  
    return py_on_sec
end

"Bringing surface of section states to statespace"
function lift(v, E, p)
    on_section(v, E, p) && return [EPS_OFF, v[1], px(v, E, p), v[2]]
end


function eqm!(du, u, p, t)
    a, m, w = p 
    x, y, px, py = u
    du[1] =  px/m
    du[2] =  py/m
    du[3] = -m*w^2*x - 2a*x*y
    du[4] = -m*w^2*y - a*(x^2 - y^2)
end



"""
Build the fast (residuals) and dense (trajectories) integrators.
 
Returns the nmax Refs as well -- always construct `SectionParams` from these
returned Refs, never from freshly-made ones, or writes to `prm.nmax_*` will be
invisible to the callbacks.
"""
function create_integrators(p; tmax = 20_000.0, dt = 0.3, nfast = 1, ndense = 40)
    nmax_fast, nmax_dense = Ref(nfast), Ref(ndense)
    
    condition(u, t, integ) = u[1]
    
    # --- fast: no trajectory saved, terminates as soon as n crossings are in ---
    yf, pyf, tsf = Float64[], Float64[], Float64[]
    affect_f!(integ) = begin
        push!(yf, integ.u[2]); push!(pyf, integ.u[4]); push!(tsf, integ.t)
        length(yf) ≥ nmax_fast[] && terminate!(integ)
    end
    integ_fast = init(ODEProblem(eqm!, zeros(4), (0.0, tmax), p),
                      Vern9(); abstol = TOL_INTEGRATE, reltol = TOL_INTEGRATE,
                      save_everystep = false, save_start = false,
                      callback = ContinuousCallback(condition, affect_f!, nothing;
                                                    abstol = TOL_CALLBACK))
 
    # --- dense: saves the trajectory for plotting / period detection ---
    yd, pyd, tsd = Float64[], Float64[], Float64[]
    affect_d!(integ) = begin
        push!(yd, integ.u[2]); push!(pyd, integ.u[4]); push!(tsd, integ.t)
        length(yd) ≥ nmax_dense[] && terminate!(integ)
    end
    integ_dense = init(ODEProblem(eqm!, zeros(4), (0.0, tmax), p),
                       Vern9(); abstol = TOL_INTEGRATE, reltol = TOL_INTEGRATE, saveat = dt,
                       callback = ContinuousCallback(condition, affect_d!, nothing;
                                                     abstol = TOL_CALLBACK))
 
    return (; integ_fast, yf, pyf, tsf, nmax_fast,
              integ_dense, yd, pyd, tsd, nmax_dense)
end

"""
The simulation parameters creats integrators, for finding periodic orbits (.integ_fast) and obtaining the periodic trajectories (.integ_dense). 
"""
function SimParams(E, p; tmax = 20_000.0, nfast = 1, ndense = 40)
    b = create_integrators(p; tmax, nfast, ndense)
    return SimParams(
        b.integ_fast, b.yf, b.pyf, b.tsf, b.nmax_fast, 
        b.integ_dense, b.yd, b.pyd, b.tsd, b.nmax_dense, 
        E, p, tmax
        )
end




"n-th crossing of the section, starting from v. This is T^n(v)."
function pMap(v::Vector{Float64}, n::Int, prm::SimParams)
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

Fres(v::Vector{Float64}, n::Int, prm::SimParams) = pMap(v, n, prm) - v

"Finite everywhere: TrustRegion probes outside the boundary and NaN poisons it."
function Fres_safe(v::Vector{Float64}, n::Int, prm::SimParams)
    a = px2(v, prm.E, prm.p)
    a <= 0 && return fill(1.0 + 100 * sqrt(-a), 2)
    return Fres(v, n, prm)
end



function jacobian!(J, v, n, prm, d)
    E, p = prm.E, prm.p
    for (j, e) in enumerate(([d, 0.0], [0.0, d]))
        vp, vm = v .+ e, v .- e
        okp, okm = on_section(vp, E, p), on_section(vm, E, p)
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



"Selfe implemented rootfinding algorithem"
function find_orbit(v0::Vector{Float64}, n::Int, prm::SimParams;
                    N_max = 100, d = 1e-7, tol = TOL_PO_ROOT, 
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


"NonlinearSolve variant. AutoFiniteDiff is mandatory: the ODE callback rejects Duals."
function solve_orbit(v0, n, prm; tol = TOL_PMAP_RES, maxiters = 300)
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


"return sol, pMap, tsd"
function section_trj(v, prm)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && error("point $v outside energy boundary")

    empty!(prm.yd); empty!(prm.pyd); empty!(prm.tsd)
    reinit!(prm.integ_dense, u0)
    solve!(prm.integ_dense)
    sol = prm.integ_dense.sol
    return sol, permutedims([prm.yd prm.pyd]), copy(prm.tsd)
end

"Find the smallest k with |T^k v - v| < tol.

returns named tuple (;traj, pMap, Nperiod, Tperiod)"
function minPeriodicity(v, prm; tol = 1e-8, search = 40)
    prm.nmax_dense[] = search          # dense integrator only — no interference with T

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


"""returns a row to the orbit

    return (; E = prm.E, n,
              prime = orb.Nperiod, T = orb.Tperiod,
              trace = τ, detDT = det(DT),
              class = KIND_LABEL[index], index = index,
              resnorm = res.resnorm,
              sec_traj = orb.pMap,
              state_traj = orb.traj,
              id=id, origin=origin)
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
 
    return (; E = prm.E, n,
              prime = orb.Nperiod, T = orb.Tperiod,
              trace = τ, detDT = det(DT),
              class = KIND_LABEL[index], index = index,
              resnorm = res.resnorm,
              sec_traj = orb.pMap,
              state_traj = orb.traj,
              id=id, origin=origin)
end





v0 = [0.2, 0.3]  # y, py
n    = 4
E    = 0.1144

prm = SimParams(E, (1.0,1.0,1.0))
orbit = analyse_seed(v0, n, prm); orbit === nothing && println("orbit failed")



# for configuration space plots
const RGRID  = range(-1.0, 1.0, length = 240)
const EPOT   = [V(x, y; p=(1.0,1.0,1.0)) for x in RGRID, y in RGRID]
const LEVELS = collect(logrange(5.0*0.009, 6.9*0.089, 7))

cfig = Figure(size=(1400,900))
cax = Axis(cfig, xlabel="x", ylabel="y")
scatter!(cax, orbit.state_traj[1:2,:])
contour!(cax, RGRID, RGRID, EPOT; LEVELS, colormap = :hsv, labels = true)

display(cfig)