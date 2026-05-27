using GLMakie

include("../models/henon_heiles.jl")
include("../analysis/poincare.jl")

# --- initial conditions ---
E  = Observable(1/8)
x0 = Observable(0.0)
y0 = Observable(0.0)

# --- derived trajectory ---
sol = lift(E, x0, y0) do e, x, y
    px = sqrt(max(0.0, 2e - y^2 - x^2 - 2x^2*y + 2y^3/3))
    solve_trajectory([x, y, px, 0.0], tspan=(0.0, 300.0))
end

traj_x  = lift(s -> [u[1] for u in s.u], sol)
traj_y  = lift(s -> [u[2] for u in s.u], sol)
ps_x, ps_px = lift(s -> poincare_section(s), sol) |> s -> (lift(x->x[1],s), lift(x->x[2],s))

# --- layout ---
fig = Figure(resolution=(1200, 700))

ax1 = Axis(fig[1,1], title="Phase space (x, y)", xlabel="x", ylabel="y")
ax2 = Axis(fig[1,2], title="Poincaré section (y=0)", xlabel="x", ylabel="pₓ")
ax3 = Axis(fig[2,1:2], title="x(t)", xlabel="t", ylabel="x")

lines!(ax1, traj_x, traj_y, linewidth=0.5)
scatter!(ax2, ps_x, ps_px, markersize=2)
lines!(ax3, lift(s->s.t, sol), traj_x, linewidth=0.5)

# --- controls ---
sg = SliderGrid(fig[3, 1:2],
    (label="Energy E",  range=0.01:0.01:0.5, startvalue=1/8),
    (label="x₀",        range=-0.4:0.01:0.4,  startvalue=0.0),
    (label="y₀",        range=-0.4:0.01:0.4,  startvalue=0.0),
)
connect!(E,  sg.sliders[1].value)
connect!(x0, sg.sliders[2].value)
connect!(y0, sg.sliders[3].value)

display(fig)
