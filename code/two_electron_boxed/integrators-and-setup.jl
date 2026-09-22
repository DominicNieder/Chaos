# This file contains functions:
#
# -----------------------------------------------------------------------------------------
# SectionParams{IF, ID, P}
# -----------------------------------------------------------------------------------------
# set_energy!(prm::SectionParams, E)
# -----------------------------------------------------------------------------------------
# get_section_p(prm::SectionParams; integ=:dense)
# -----------------------------------------------------------------------------------------
#get_section_t(prm::SectionParams; integ=:dense)
# -----------------------------------------------------------------------------------------
# eom!()
# -----------------------------------------------------------------------------------------
# wall_cb(p, s, section, pts; cc_tol = CC_TOL)
# -----------------------------------------------------------------------------------------
# wall_callback(p; section = (1, 1), cc_tol=CC_TOL)
# -----------------------------------------------------------------------------------------
# integrator(p; 
#            tmax = 20_000.0, kind=:fast, n = 1, cc_tol = CC_TOL, 
#            int_tol=INT_TOL, save_everystep = false, save_start = false)
# -----------------------------------------------------------------------------------------
# create_integrators(p; 
#              tmax = 20_000.0, nfast = 1, ndense = 40, 
#              cc_tol = CC_TOL, int_tol = INT_TOL, 
#              save_everystep = save_everystep, save_start = save_start)
# -----------------------------------------------------------------------------------------
# SectionParams(E, p, num_int; 
#           tmax = 20_000.0, nfast = 1, ndense=40, 
#           cc_tol = CC_TOL, int_tol = INT_TOL, save_everystep = true, save_start = true)
# -----------------------------------------------------------------------------------------

mutable struct SectionParams{IF, ID, P}
    integ_fast  :: IF
    ptsf        :: Vector{Tuple{Float64, Float64, Float64}}
    nmax_fast   :: Ref{Int}
 
    integ_dense :: ID
    ptsd        :: Vector{Tuple{Float64, Float64, Float64}}
    nmax_dense  :: Ref{Int}

    E           :: Float64
    p           :: P
    tmax        :: Float64
end

set_energy!(prm::SectionParams, E) = (prm.E = E) 

get_section_p(prm::SectionParams; integ=:dense)  = integ==:dense ? get_section_p(prm.ptsd) : get_section_p(prm.ptsf)
get_section_p(pts::Vector{Tuple{Float64, Float64, Float64}}) = first.(pts,2)

get_section_t(prm::SectionParams; integ=:dense)  = integ==:dense ? get_section_t(prm.ptsd) : get_section_t(prm.ptsf)
get_section_t(pts::Vector{Tuple{Float64, Float64, Float64}})  = last.(pts)


"""
    Equations of motion for two-boxed-electron model
"""
function eom!(du, u, p, t)
    x1, x2, p1, p2 = u
    d = x1 - x2
    f = p.C*sign(d) / d^2        # force on particle 1
    du[1] = p1 / p.m1
    du[2] = p2 / p.m2
    du[3] =  f
    du[4] = -f 
end


"""

    pts: [(sec_y, sec_py, sec_t),...]
    wall_callback(p; section = (1, 2))
 
Elastic reflection at all four walls. Every time the wall given by `section`
is hit, the *other* particle's (x, p) is recorded. Returns `(cb, pts)`.
"""
function wall_callback(p; cc_tol = CC_TOL)
    pts     = Tuple{Float64,Float64,Float64}[]
    boxes   = get_boxes(p)
    section = p.C > 0 ? (1,1) : (1,2)
    walls   = [(1,1), (1,2), (2,1), (2,2)]   # event index -> (particle i, side s)

    function condition!(out, u, t, integrator)
        for (k, (i, s)) in enumerate(walls)
            out[k] = u[i] - boxes[i][s]
        end
    end

    function bounce!(integrator, k::Int)
        i, s    = walls[k]
        outward = s == 1 ? -1.0 : 1.0
        mom_idx = 2 + i
        if sign(integrator.u[mom_idx]) == outward
            integrator.u[mom_idx] = -integrator.u[mom_idx]
            if (i, s) == section
                j = 3 - i
                push!(pts, (integrator.u[j], integrator.u[2+j], integrator.t))
            end
        end
    end

    # k can be a plain Int (single event) or a sign/flag vector (tied events at the
    # same instant) -- handle both.
    function affect!(integrator, k)
        if k isa Integer
            bounce!(integrator, k)
        else
            for idx in eachindex(k)
                k[idx] != 0 && bounce!(integrator, idx)
            end
        end
    end

    cb = VectorContinuousCallback(condition!, affect!, 4; interp_points = 20, abstol = cc_tol)
    return cb, pts
end



"""
    p=(;C, m1, m2, del)

    return (; integ, pts)

pts := surface of section

Build the fast (residuals) and dense (trajectories) integrators.
 
Returns the nmax Refs as well -- always construct `SectionParams` from these
returned Refs, never from freshly-made ones, or writes to `prm.nmax_*` will be
invisible to the callbacks.
"""
function integrator(p; 
    tmax = 20_000.0, kind=:fast, n = 1, cc_tol = CC_TOL, int_tol=INT_TOL, save_everystep = false, save_start = false)
    nmax = Ref(n)
 
    # --- fast: no trajectory saved, terminates as soon as n crossings are in ---

    callback, pts   = wall_callback(p; cc_tol=cc_tol)
    prob            = ODEProblem(eom!, zeros(4), (0.0, tmax), p)
    if kind == :dense
        integ      = init(prob, Vern9(); 
                            abstol = int_tol, reltol = int_tol,
                            save_everystep = save_everystep, save_start = save_start,
                            callback = callback)
    elseif kind == :fast
        integ      = init(prob, Vern9();
                            abstol = int_tol, reltol = int_tol,
                            save_everystep = false, save_start = false,
                            callback = callback)
    else
        error("kind must be :fast or :dense")
    end
    return (; integ, pts, nmax)
end


"""
    retrun (; integ_fast = i_fast.integ, secf=first.(i_fast.pts,2), tsf = last.(i_fast.pts), nmax_fast = i_fast.nmax, 
    integ_dense = i_dense.integ, secd = first.(i_dense.pts,2), tsd = last.(i_dense.pts), nmax_dense = i_dense.nmax)

    pts = Tuple{Float64, Float64, Float64}[(y, py, t),...]

Build the fast (residuals) and dense (trajectories) integrators.
 
If you want to look at the orbits in configuration space, set save_everystep = true in the dense integrator, and then use `section_trj` to get the trajectory.
"""
function create_integrators(p; 
             tmax = 20_000.0, nfast = 1, ndense = 40, 
             cc_tol = CC_TOL, int_tol = INT_TOL, 
             save_everystep = save_everystep, save_start = save_start)

    # --- fast: no trajectory saved, terminates as soon as n crossings are in ---
    i_fast = integrator(p; tmax, kind=:fast, n=nfast, cc_tol, int_tol, save_everystep = false, save_start = false)
    # --- dense: saves the trajectory for plotting / period detection ---
    i_dense = integrator(p; tmax, kind=:dense, n=ndense, cc_tol, int_tol, save_everystep = save_everystep, save_start = save_start)
    return (; integ_fast = i_fast.integ, ptsf = i_fast.pts, nmax_fast = i_fast.nmax, 
    integ_dense = i_dense.integ, ptsd = i_dense.pts, nmax_dense = i_dense.nmax)
end


function SectionParams(E, p; 
          tmax = 20_000.0, nfast = 1, ndense=40, 
          cc_tol = CC_TOL, int_tol = INT_TOL, save_everystep = true, save_start = true)

    b = create_integrators(p; tmax, nfast=nfast, ndense=ndense, cc_tol, int_tol, save_everystep = save_everystep, save_start = save_start)
    return SectionParams(b.integ_fast, b.ptsf, b.nmax_fast, b.integ_dense, b.ptsd, b.nmax_dense, E, p, tmax)
end