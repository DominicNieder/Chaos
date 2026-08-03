module trajectories

"""
makes a function that filters the surface of section values whilst solving the ODE
"""
function section_callback(section_q, section_p)
    condition(u, t, integrator) = u[1]
    affect!(integrator) = begin
        if integrator.u[3] >= 0
            push!(section_q, integrator.u[2])
            push!(section_p, integrator.u[4])
        end
    end
    ContinuousCallback(condition, affect!, nothing)
end


"""
u0 -> initial conditions
tspan   =   (tmin, tmax)
dt      =   saved time points (adaptive integration)
p       =   parameters of the henon heiles modle

Uses Vern9 algorithem: 
"Verner's “Most Efficient” 9/8 Runge-Kutta method. (lazy 9th order interpolant)"
"""
function solve_trajectory(
    u0,
    p, 
    tspan, 
    dt; 
    callback=nothing
    )
    prob = ODEProblem(
        equations!, u0, tspan, p
        )
    solve(
        prob, 
        Vern9(), 
        abstol=1e-9, 
        reltol=1e-9, 
        saveat=dt, 
        callback=callback
        )
end




end # module