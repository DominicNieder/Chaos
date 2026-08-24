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
const FIG_DIR    = joinpath(@__DIR__, "../../figures/henon-heiles/periodic-orbits/")
data_file        = joinpath(DATA_DIR, "E0.1127-T10000.0-py0.0-n256.jld2")

configurations = JSON3.read(read(CONFIG_DIR, String))
cfg            = configurations.explore
param          = [Float64(cfg.a.value), Float64(cfg.m.value), Float64(cfg.w.value)]
dt             = Float64(cfg.dt.value)
x0             = Float64(cfg.x0.value)

struct SectionParams{IF, ID, P}
    integ_fast  :: IF
    yf :: Vector{Float64};  pyf :: Vector{Float64};  tsf :: Vector{Float64}
    integ_dense :: ID
    yd :: Vector{Float64};  pyd :: Vector{Float64};  tsd :: Vector{Float64}
    nmax_fast  :: Base.RefValue{Int}
    nmax_dense :: Base.RefValue{Int}
    E    :: Float64
    p    :: P
    tmax :: Float64
end

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
    id::Int
end

# --- background section ---
data = load(data_file, "results")
y_all, py_all = Float32[], Float32[]
for d in data
    append!(y_all,  d.sec_y)
    append!(py_all, d.sec_py)
end

const EPS_OFF = -1e-15

"""
returns -> 1 := elliptic, 2 := hyperbollic (stable), 3 := hyperbollic (unstable), 4 := parabolic
"""
function kind_index(τ; ε = 1e-6)
    abs(abs(τ) - 2) < ε && return 4   # parabolic
    abs(τ) < 2 && return 1            # elliptic
    τ > 0 ? 2 : 3                     # hyperbollic (unstable, stable)
end

const KIND_LS    = [:solid, :dash, :dashdot, :dot]
const KIND_MS    = [:circle, :xcross, :diamond, :utriangle]
const KIND_LABEL = ["elliptic", "hyperbolic", "inverse hyperbolic", "parabolic"]


px2(y, py, E, p)    = 2p[2] * (E - HenonHeiles.potential(0.0, y, p)) - py^2
in_section(v, E, p) = px2(v[1], v[2], E, p) > 0



"Lift (y, py) to a 4D state. Offset follows sign(px) -> no phantom t=0 crossing."
function lift(v, E, p; sgn = +1)
    a = px2(v[1], v[2], E, p)
    a <= 0 && return nothing
    return [EPS_OFF, v[1], sgn * sqrt(a), v[2]]
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

"""
    Find the smallest k with |T^k v - v| < tol.

returns named tuple (;traj, pMap, Nperiod, Tperiod)
"""
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

"""
    Two integrators, one for finding periodic orbit (integ_fast) and another (integ_dense) of for finding minimal periodicity and finding the trajectories of the orbit.

return (; integ_fast, yf, pyf, tsf, integ_dense, yd, pyd, tsd)
"""
function create_integrators(p; tmax = 20_000.0, nmax = Ref(1))
    # fast: for residual evaluations, no trajectory saved
    yf, pyf, tsf = Float64[], Float64[], Float64[]
    condition(u, t, integ) = u[1]
    affect!(integ) = begin
        push!(yf, integ.u[2]); push!(pyf, integ.u[4]); push!(tsf, integ.t)
        length(yf) ≥ nmax[] && terminate!(integ)
    end
    cbf = ContinuousCallback(condition, affect!, nothing; abstol = 1e-13)
    integ_fast = init(ODEProblem(HenonHeiles.equations!, zeros(4), (0.0, tmax), p),
                      Vern9(); abstol = 1e-14, reltol = 1e-14,
                      save_everystep = false, save_start = false, callback = cbf)

    # dense: for plotting
    yd, pyd, tsd = Float64[], Float64[], Float64[]
    affect2!(integ) = begin
        push!(yd, integ.u[2]); push!(pyd, integ.u[4]); push!(tsd, integ.t)
        length(yd) ≥ nmax[] && terminate!(integ)
    end
    cbd = ContinuousCallback(condition, affect2!, nothing; abstol = 1e-13)
    integ_dense = init(ODEProblem(HenonHeiles.equations!, zeros(4), (0.0, tmax), p),
                       Vern9(); abstol = 1e-14, reltol = 1e-14,
                       saveat = dt, callback = cbd)

    return (; integ_fast, yf, pyf, tsf, integ_dense, yd, pyd, tsd)
end


""
function T(v, n::Int, prm::SectionParams)
    u0 = lift(v, prm.E, prm.p)
    u0 === nothing && @warn "point $v outside energy boundary"
    
    n > prm.nmax_fast[] && (prm.nmax_fast[] = n)     # raise the cap if this n needs more

    empty!(prm.yf); empty!(prm.pyf); empty!(prm.tsf)
    reinit!(prm.integ_fast, u0)
    solve!(prm.integ_fast)
    length(prm.yf) < n && error("only $(length(prm.yf)) crossings (need $n)")
    return [prm.yf[n], prm.pyf[n]]
end

Fres(v, n, prm) = T(v, n, prm) - v

function Fres_safe(v, n, prm)
    a = px2(v[1], v[2], prm.E, prm.p)
    a <= 0 && return fill(1.0 + 100*sqrt(-a), 2)   # penalty scales with violation
    return Fres(v, n, prm)
end

function jacobian!(J, v, n, prm, d)
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
            @warn ("both probes outside boundary at $v — step too large or v on the edge")
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
                    N_max = 100, d = 1e-7, tol = 1e-11, 
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

function solve_orbit(v0, prm; tol = 1e-11, maxiters = 300)
    v = collect(float.(v0))
    in_section(v, prm.E, prm.p) ||
        return (v = v, DT = nothing, converged = false,
                resnorm = Inf, history = zeros(2,0), comment = "seed outside boundary")

    prob = NonlinearProblem((w, q) -> Fres_safe(w, q.n, q), v, prm)    
    sol  = solve(prob, TrustRegion(autodiff = AutoFiniteDiff(fdtype = Val(:central)));
                 abstol = tol, maxiters)

    ok = SciMLBase.successful_retcode(sol) && in_section(sol.u, prm.E, prm.p)
    ok || return (v = sol.u, DT = nothing, converged = false,
                  resnorm = norm(sol.resid), history = zeros(2,0),
                  comment = "$(sol.retcode)")

    DT = jacobian!(zeros(2,2), sol.u, n, prm, 1e-7) + I
    return (v = sol.u, DT = DT, converged = true,
            resnorm = norm(sol.resid), history = zeros(2,0),
            comment = "$(sol.retcode)  det(DT) = $(det(DT))")
end

"Is the orbit invariant under S(y,py) = (y,-py)?"
function is_sym(sec; tol = 1e-6)
    Ssec = sec .* [1.0, -1.0]
    for j in axes(Ssec, 2)
        any(i -> norm(Ssec[:, j] - sec[:, i]) < tol, axes(sec, 2)) || return false
    end
    return true
end

P_py(v) = [v[1], -v[2]]
P_pypx(u) = [u[1], u[2], -u[3], -u[4]]
P_x(u) = [-u[1], u[2], -u[3], u[4]]


"rotate momenta and position by an angle θ"
function rotate_space(v, θ)
    c, s = cos(θ) ,sin(θ)
    return [c -s; s c]*v 
end

function rotate_state(u, θ)
    c, s = cos(θ) ,sin(θ)
    return [c*u[1]-s*u[2], s*u[1]+c*u[2], c*u[3]-s*u[4], s*u[3]+c*u[4]]
end

function is_rot_sym(sec, E, p; tol = 1e-6)
    Ssec = rotate_state(lift(sec, E, p))
    for j in axes(Ssec, 2)
        any(i -> norm(Ssec[:, j] - sec[:, i]) < tol, axes(sec, 2)) || return false
    end
    return true
end

function rotated_root(v, prm, θ)
    u_rot = rotate_state(lift(v, prm.E, prm.p), θ)
    empty!(prm.yd); empty!(prm.pyd); empty!(prm.tsd)
    reinit!(prm.integ_dense, u_rot)
    solve!(prm.integ_dense)
    isempty(prm.yd) && error("rotated state never crossed the section")
    return [prm.yd[1], prm.pyd[1]]
end

" 
E::Float64; n::Int
seed_y::Float64; seed_py::Float64
y::Float64; py::Float64
prime::Int; T::Float64
trace::Float64; detDT::Float64
resnorm::Float64; iters::Int
sec_y::Vector{Float64}; sec_py::Vector{Float64}
history_y::Vector{Float64}; history_py::Vector{Float64}
traj_x::Vector{Float64}; traj_y::Vector{Float64}
"
orbit_table() = DataFrame(Row[])

"returns a row to the orbit"
function analyse_seed(v0, n, prm; id=0)::Union{Row,Nothing}
    res = find_orbit(v0, n, prm)
    res.converged || return nothing
    orb = minPeriodicity(res.v, prm; tol = 1e-8)
    orb.Nperiod === nothing && return nothing
    DT=get_DT(res.v, orb.Nperiod, prm)
    trace=tr(DT); detDT=det(DT)
    index =  kind_index(trace; ε = 1e-6)
    return (; E = prm.E, n = n, seed_y = v0[1], seed_py = v0[2],
                     y = res.v[1], py = res.v[2],
                     prime = orb.Nperiod, T = orb.Tperiod,
                     trace = trace, detDT = detDT,
                     resnorm = res.resnorm, iters = size(res.history, 2),
                     sec_y = orb.pMap[1,:], sec_py=orb.pMap[2,:],
                     history_y = res.history[1,:], history_py = res.history[2,:],
                     traj_x = [u[1] for u in orb.traj], traj_y = [u[2] for u in orb.traj],
                     class = KIND_LABEL[index], index=index, id=id)
end



p       =   (1.0, 1.0, 1.0)
E0      =   0.1127
nmax_search = 40
nmax        =   1
tmax    =   20_000.0
y0, py0, c0 =    0.0, 0.2, C_TEAL
y1, py1, c1 =   -0.35, 0.0, C_ORANGE
y2, py2, c2 =    0.27, 0.3, C_GREEN
y3, py3, c3 =    0.3, 0.0, C_PURPLE
v_init      =   [[y0, py0], [y1, py1], [y2, py2], [y3, py3]] #, 
color_scheme= [c0, c1, c2, c3]

bff = create_integrators(p; nmax=Ref(nmax_search))

prm = SectionParams(bff.integ_fast, bff.yf, bff.pyf, bff.tsf,
                         bff.integ_dense, bff.yd, bff.pyd, bff.tsd,
                         Ref(nmax), Ref(nmax_search), E0, p, tmax)

orbits = orbit_table()
for (i,v0) in enumerate(v_init)
    orbit = analyse_seed(v0, nmax, prm; id=i)
    orbit === nothing || push!(orbits, orbit)
end
# now I want to see how the symmetry transformations effect my orbits

tasks= ["Rotation_sym1", "Rotation_sym2", "Rotation_sym3"]



base = collect(eachrow(orbits))
# for i in 1:2, orb in base
#     v1 = try
#         rotated_root([orb.y, orb.py], prm, i*2π/3)
#     catch e
#         @warn "rotation failed" seed=(orb.y, orb.py) exception=e
#         continue
#     end
#     @show i*2/3, norm(Fres(v1, orb.prime, prm)), orb.id
#     new = analyse_seed(v1, 1, prm; id= orb.id)
#     # new.id=orb.id
#     new === nothing || push!(orbits, new)
# end

for orb in base
    sec = permutedims([orb.sec_y orb.sec_py])     # rebuild the 2 x k matrix
    is_sym(sec; tol = 1e-6) && continue            # already closed under S — skip
    v1  = P_py([orb.y, orb.py])
    new = analyse_seed(v1, orb.prime, prm; id = orb.id)
    new === nothing || push!(orbits, new)
end




r      = range(-1.0, 1.0, length=120)
levels = logrange(5.0*0.009, 6.9*0.089, 7)
epot   = [HenonHeiles.potential(x,y, param) for x in r, y in r]

fig_conf = Figure(size = (1400, 900))
ax_conf  = Axis(fig_conf[1, 1], xlabel = "x", ylabel = "y",
             title = "config space, n = 1", aspect = DataAspect())
contour!(ax_conf, r, r, epot, labels=true, levels=levels, colormap=:hsv, colorscale=identity)
for (j,o) in enumerate(eachrow(orbits))
    # j >4 && lines!(ax_conf, o.traj_x, o.traj_y, color = color_scheme[o.id],
    #        label = "T=$(round(o.T; digits=1)), $(o.class)", linestyle=KIND_LS[o.index])
    # j==4 && lines!(ax_conf, o.traj_x, o.traj_y, color = color_scheme[o.id],
    #        label = "T=$(round(o.T; digits=1)), $(o.class)", linestyle=(:dash))
    lines!(ax_conf, o.traj_x, o.traj_y, color = color_scheme[o.id],
            label = "T=$(round(o.T; digits=1)), $(o.class)", linestyle=KIND_LS[o.index])
end #round(x; digits = 3)
Legend(fig_conf[2, 1], ax_conf; orientation = :horizontal, framevisible = false)
display(GLMakie.Screen(), fig_conf)

fig_p = Figure(size=(1400,900))
ax_p = Axis(fig_p[1, 1], xlabel = "y", ylabel = "p_y",
               title = "section, n = 1, prime = 1")
    

y_max, py_max = HenonHeiles.section_boundary_ranges(E0, param, 120)
boundary      = HenonHeiles.section_boundary(y_max, py_max) 
scatter!(ax_p, y_all, py_all, color = (:grey, 0.5), markersize = 1.5)
scatter!(ax_p, boundary, color = C_CREAM, markersize = 4)
# scatterlines!(ax_p, res.history[1, :], res.history[2, :],
                    # color = C_GOLD, markersize = 8, label = "Newton path")
j=1
for (j, o) in enumerate(eachrow(orbits))
    # j>4 && scatter!(ax_p, o.sec_y, o.sec_py,
    #          color = color_scheme[o.id], markersize = 12,
    #          marker = KIND_MS[o.index],
    #          label = KIND_LABEL[o.index])
    scatter!(ax_p, o.sec_y, o.sec_py,
             color = color_scheme[o.id], markersize = 12,
             marker = KIND_MS[o.index],
             label = KIND_LABEL[o.index])
end
Legend(fig_p[2, 1], ax_p; orientation = :horizontal, framevisible = false)
display(GLMakie.Screen(), fig_p)




# sol=Orbit_finder(v0, n, E0; theta=2/3)

save(joinpath(FIG_DIR, "phsp-two-orbits-by-Sreflection-AndBase-$(E0)-n$nmax.png"),fig_p; px_per_unit = 2)
save(joinpath(FIG_DIR, "conf-two-orbits-by-Sreflection-AndBase-$(E0)-n$nmax.png"), fig_conf; px_per_unit = 2)

println("\n")
for (j,o) in enumerate(eachrow(orbits))

    j <= 4 && println("$j) Base orbit: $(o.id)\n y=$(o.y), py=$(o.py)\n|r|=$(o.resnorm)")
end
