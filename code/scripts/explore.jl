using GLMakie

include("../models/henon_heiles.jl")
include("../analysis/poincare.jl")

# --- initial conditions ---
E  = Observable(init_E0)
x0 = Observable(init_x0)
y0 = Observable(init_y0)
py0= Observable(init_py0)

py0_eff = lift(E, x0, y0, py0) do e, x, y, py
    a, m, w = param
    v    = m*w^2/2*(x^2 + y^2) + a*(x^2*y - y^3/3)
    pmax = sqrt(max(0.0, 2*e*m - v))
    clamp(py, 0.0, pmax)
end

# --- derived trajectory ---
sol = lift(E, x0, y0, py0_eff) do e, x, y, py
    a, m, w = param
    diff = 2e*m - m*w^2/2*(y^2 + x^2) - a*(x^2*y - y^3/3) - py^2
    px   = sqrt(max(0.0, diff))
    solve_trajectory([x, y, px, py], tspan=(0.0, num_int[1]), dt=num_int[2], p=param)
end

traj_x  = lift(s -> [u[1] for u in s.u], sol)
traj_y  = lift(s -> [u[2] for u in s.u], sol)
ps_x, ps_px = lift(s -> poincare_section(s), sol) |> s -> (lift(x->x[1],s), lift(x->x[2],s))

# --- layout ---
fig = Figure(resolution=(1200, 700))

ax1 = Axis(fig[1,1], title="Space (x, y)", xlabel="x", ylabel="y")
ax2 = Axis(fig[1,2], title="Poincaré section (y=0)", xlabel="x", ylabel="pₓ")
ax3 = Axis(fig[2,1:2], title="x(t)", xlabel="t", ylabel="x")

lines!(ax1, traj_x, traj_y, linewidth=0.5)
scatter!(ax2, ps_x, ps_px, markersize=2)
lines!(ax3, lift(s->s.t, sol), traj_x, linewidth=0.5)

# --- controls ---
sg = SliderGrid(fig[3, 1:2],
    (label="Energy E",  range=0.01:0.000001:0.9,startvalue=0.08333),
    (label="x₀",        range=-0.4:0.01:0.4,    startvalue=0.0),
    (label="y₀",        range=-0.4:0.01:0.4,    startvalue=0.0),
    (label="py₀",       range=0.0:0.01:0.4,     startvalue=0.0)
)
connect!(E,     sg.sliders[1].value)
connect!(x0,    sg.sliders[2].value)
connect!(y0,    sg.sliders[3].value)
connect!(py0,   sg.sliders[4].value)

Label(fig[4, 1:2], lift(v -> "effective py₀ = $(round(v, digits=4))", py0_eff))

display(fig)
