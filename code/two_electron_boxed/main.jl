include("two_electron_boxed_apply_monodromy.jl")


BLAS.set_num_threads(1)          # your linear algebra is 2x2; BLAS threads only compete


# ============================================================================================
#                       Save and data directories
# ============================================================================================

const DATA_DIR   = joinpath(@__DIR__, "../../data/two-boxed-charges/")
const FIG_DIR    = joinpath(@__DIR__, "../../figures/two-boxed-charges/")



# ============================================================================================
#                       constants for precission of the simulations
# ============================================================================================


const CC_TOL  = 1e-13
const INT_TOL = 1e-14          # how precise the integrator should be for sattisfactory convergance
const PMAP_ROOT_TOL = 1e-11    # poincare returnmap precission to find periodic orbit, i.e. the root
const PMAP_PRIME_TOL = 1e-9    # tollerace for saying that two roots belong to the same orbit
const DEL_BOX = 1e-8           # offset between the two boxes that breaks the symmetry of the charged point particles


p = (;C=-1, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= DEL_BOX)
E  = 500 / (p.C/p.L)

scatter(boundary(E;p=p,n=10000))
seeds = [[0.0006,0.0], [0.0006,20]]  # C<0
#seeds =  [[0.2,0.0], [0.32,0.0], [0.5,0.0]]
u0s = [lift(v,E,p) for v in seeds]



int_time = 10_000



set_style!(:dark)

fig = Figure(size=(1300,900))
ax_sec    = Axis(fig[2,1:2], xlabel=L"x_2\, [L]", ylabel=L"p_2")
ax_traj   = Axis(fig[1,1:2], xlabel=L"\hat{x}_1\, [\text{L}]", ylabel=L"\hat{x}_2 \, [\text{L}]")
ax_energy = Axis(fig[2,3],   xlabel=L"t", ylabel=L"\text{Energy}\,[\mathcal{C}/L]")
# [hlines!(ax_traj, get_boxes(p)[j][i], color=INK, linestyle=:dash) for i in 1:2, j in 1:2]

Es0 = Float64[]                                  # initial energy per seed, for the legend
for (i, u0) in enumerate(u0s)
    E0 = energy(u0, p); push!(Es0, E0)
    u, t, pts = get_traj(u0, int_time; p=p, abstol=INT_TOL, reltol=INT_TOL)

    xs   = first.(u, 2)
    xone = Point2f.(zip(t, first.(xs)))
    xtwo = Point2f.(zip(t, last.(xs)))
    sec  = Point2f.(first.(pts, 2))

    lines!(ax_traj, Point2f.(xs), color=(pick_color(i),0.8))
    # lines!(ax_traj, xtwo[1:2000], color=pick_color(i), linestyle=:dash)
    scatter!(ax_sec, sec, color=pick_color(i), markersize=4)
    lines!(ax_energy, t, [abs(energy(ui,p) - E)*p.L/E/p.C for ui in u], color=pick_color(i))
end
scatter!(ax_sec, boundary(E; p,n=1000), color=INK, markersize=3)

# ── legend: component (line style) + seeds (colour, with E₀) ──
comp_elems = [LineElement(color=:gray70, linestyle=:solid, linewidth=2),
              LineElement(color=:gray70, linestyle=:dash,  linewidth=2),
              LineElement(color=INK,     linestyle=:dash,  linewidth=1.5),
              MarkerElement(color=INK, marker=:circle, markersize=6)]
comp_labels = [L"x_1", L"x_2", "box edges", "section boundary"]

seed_elems  = [LineElement(color=pick_color(i), linewidth=3) for i in eachindex(u0s)]
seed_labels = ["seed $(seeds[i])   " for i in eachindex(u0s)]

Legend(fig[1,3], [comp_elems, seed_elems], [comp_labels, seed_labels],
       ["Component", "Seeds  (E = $(round(E, sigdigits=4)))"];
       framevisible   = false,
       tellheight     = false,
       halign         = :left, valign = :top,
       titlehalign    = :left, titlesize = 12,
       gridshalign    = :left,
       labelsize      = 11,
       patchsize      = (18, 10),
       rowgap         = 1,
       groupgap       = 8,
       padding        = (2, 2, 2, 2))

colsize!(fig.layout, 3, Auto(0.45))   # keep column 3 narrow so the legend sits snug

folder = joinpath(FIG_DIR, "explore-surface-sections/")
#save(joinpath(folder, "some_traj@E$E.png"), fig)   # save AFTER the legend exists
display(fig)