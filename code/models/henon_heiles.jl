module HenonHeiles

using OrdinaryDiffEq, JSON3, Polynomials


export load_config, save_config, equations!, solve_trajectory, potential, kinetic,
       hamiltonian, section_callback, limit_of_initial_y0, limit_of_initial_py0, section_boundary_ranges, section_boundary, SimConfig


struct SimConfig
    a::Float64; m::Float64; w::Float64
    E::Float64; x0::Float64; y0::Float64; py0::Float64
    T::Float64; dt::Float64
end

function load_config(path)
    cfig = JSON3.read(read(path, String))
    cfg     = cfig.explore
    SimConfig(cfg.a.value, cfg.m.value, cfg.w.value,
              cfg.E.value, cfg.x0.value, cfg.y0.value, cfg.py0.value,
              cfg.T.value, cfg.dt.value)
end


function save_config(path, new_key, cfg::SimConfig; comment="")
    orient_file = isfile(path) ? JSON3.read(read(path, String), Dict{String,Any}) : Dict{String,Any}()

    if haskey(orient_file, new_key)
        return save_config(path, "_rename_"*new_key, cfg, comment=comment)
    else
        cfg_dict = Dict(
                "a"         =>  cfg.a,
                "m"         =>  cfg.m,
                "w"         =>  cfg.w,
                "T"         =>  cfg.T,
                "dt"        =>  cfg.dt,
                "E"         =>  cfg.E,
                "x0"        =>  cfg.x0,
                "y0"        =>  cfg.y0,
                "py0"       =>  cfg.py0,
                "comment"   =>  comment
        )
        get!(orient_file, new_key, cfg_dict)
    end
    open(path, "w") do io
        JSON3.pretty(io, orient_file)
    end
end

"""
Equation of motion for Hénon-Heiles potential. 
"""
function equations!(du, u, p, t)
    a, m, w = p 
    x, y, px, py = u
    du[1] =  px/m
    du[2] =  py/m
    du[3] = -m*w^2*x - 2a*x*y
    du[4] = -m*w^2*y - a*(x^2 - y^2)
end


# Energies, param= (a, m, w)
potential(x, y, p) = p[2]*p[3]^2/2 *(x^2+y^2) + p[1]*(x^2*y - y^3/3)

kinetic(px, py, p) = (px^2 + py^2)/(2*p[2])

hamiltonian(x, y, px, py, p) = kinetic(px, py, p) + potential(x,y, p)


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


"""
calculates the limit for y0, for surface section x=0.
    Beware that the expression is set at
    alpha   =   1
    m       =   1
    w       =   1
"""
limit_of_initial_y0(E::Float64, p) = roots(
    Polynomial([3E/p[1], 0,-3*p[2]*p[3]^2/(2p[1]), 1])
    )


function limit_of_initial_py0(y::Real, E::Float64, p)
   sqrt(max(0.0, 2*(E - potential(0.0, y, p))))
end

function limit_of_initial_py0(y::AbstractVector{<:Real}, E::Float64, p)
    map(yi -> limit_of_initial_py0(yi, E, p), y)
end

function section_boundary_ranges(e,  param, resolution)
    y0_roots        = real.(filter(r -> abs(imag(r)) < 1e-10, limit_of_initial_y0(e, param)))
    y0_1, y0_2 = y0_roots[1], y0_roots[2]

    y0_range  = collect(range(y0_1, y0_2, length=resolution))
    py0_range = limit_of_initial_py0(y0_range, e, param)
    y0_range, py0_range
end

function section_boundary(y0_range, py0_range)
    boundary = vcat(
        collect(zip(y0_range, py0_range)),
        collect(zip(reverse(y0_range), -reverse(py0_range)))
    )
    push!(boundary, boundary[1])  # close the ring
    boundary
end

end  # module