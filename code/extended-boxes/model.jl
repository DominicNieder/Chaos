"""
    Equations of motion for two-boxed-electron model
"""
function eom!(du, u, p, t)
    x1, x2, x3, p1, p2, p3 = u
    d12 = x1 - x2
    d23 = x2 - x3
    d13 = x1 - x3
    f12 = p.e1*p.e2*sign(d12) / d12^2        # force on particle 1
    f23 = p.e2*p.e3*sign(d23) / d23^2
    f13 = p.e1*p.e3*sign(d13) / d13^2
    du[1] = p1 / p.m1
    du[2] = p2 / p.m2
    du[3] = p3 / p.m3
    du[4] =  f12 + f13
    du[5] = -f12 + f23
    du[6] = -f23 - f13
end

function dB(p)
    return ((-1.5,-0.5),(-0.5,0.5),(0.5,1.5))
end

function wall_callback(p; cc_tol = CC_TOL)
    pts     = Tuple{Float64,Float64,Float64}[]
    boxes   = boundary_conditions(p)
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



prob = prob()