using OrdinaryDiffEq, LinearAlgebra

include("../models/henon_heiles.jl")

# Variational equations for the tangent vector alongside the trajectory
function variational!(du, u, p, t)
    # u = [x, y, px, py, dx, dy, dpx, dpy]
    x, y, px, py   = u[1:4]
    dx, dy, dpx, dpy = u[5:8]

    # equations of motion
    du[1] = px;  du[2] = py
    du[3] = -x - 2x*y
    du[4] = -y - x^2 + y^2

    # Jacobian * tangent vector
    J = [ 0  0  1  0;
          0  0  0  1;
         -1-2y  -2x  0  0;
         -2x   -1+2y  0  0]
    dv = J * [dx, dy, dpx, dpy]
    du[5:8] .= dv
end

function max_lyapunov(u0; tspan=(0.0, 1000.0), dt=1.0)
    d0 = [1.0, 0.0, 0.0, 0.0]
    u_ext = vcat(u0, d0 ./ norm(d0))
    prob = ODEProblem(variational!, u_ext, tspan)
    sol = solve(prob, Vern9(), abstol=1e-10, reltol=1e-10, saveat=dt)

    λ_sum = 0.0
    for i in 2:length(sol.t)
        d = sol[5:8, i]
        λ_sum += log(norm(d))
        sol.u[i][5:8] ./= norm(d)
    end
    λ_sum / tspan[2]
end
