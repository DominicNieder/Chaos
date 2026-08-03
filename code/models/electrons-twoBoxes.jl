using OrdinaryDiffEq, JSON3, Polynomials

"""
Equation of motion for two one paticle boxs with EM interaction
"""
function equations1d!(du,u,p,t)
    C, m =p
    x1, x2, p1, p2 = u
    f = C * sign(x1 - x2) / (x1 - x2)^2
    du[1] = p1/m
    du[2] = p2/m
    du[3] = f
    du[4] = -f
end




potential2D(x1, y1, x2, y2, p) = p[1]/sqrt((x1-x2)^2 + (y1-y2)^2)

potential1D(x1,x2, p) = p[1]/sqrt((x1-x2)^2)

kinetic1D(p,param) = p^2/(2p[2])

kinetic(px1, py1, p) =kinetic1D(px1,p) kinetic1D(px2,p)

hamiltonian(x, y, px, py, p) = kinetic(px, py, p) + potential(x,y, p)

"""
Callback function that reflects at wall collisions
"""
function wall_callback(boxes, section)   # section = (i_wall, s_wall) to record on
    walls = [(i, s) for i in 1:2 for s in 1:2]
    pts = Vector{Tuple{Float64,Float64}}()   # collected section points
    function condition(out, u, t, integrator)
        for (k, (i, s)) in enumerate(walls)
            out[k] = u[i] - boxes[i][s]
        end
    end
    function affect!(integrator, k)
        i, s = walls[k]
        outward = (s == 1) ? -1 : +1
        if sign(integrator.u[2 + i]) == outward
            integrator.u[2 + i] *= -1
            if (i, s) == section
                j = 3 - i                          # the other particle
                push!(pts, (integrator.u[j], integrator.u[2 + j]))
            end
        end
    end
    cb = VectorContinuousCallback(condition, affect!, length(walls))
    return cb, pts
end

cb, pts = wall_callback(boxes, (1, 2))   # section: particle 1 hits right wall
sol = solve(prob, Vern9(); callback = cb, abstol=1e-9, reltol=1e-9)
scatter(first.(pts), last.(pts))





