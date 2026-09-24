include("two_electron_boxed_apply_monodromy.jl")


BLAS.set_num_threads(1)          # your linear algebra is 2x2; BLAS threads only compete


# ============================================================================================
#                       Save and data directories
# ============================================================================================

const CONFIG_DIR = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
const DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
const FIG_DIR    = joinpath(@__DIR__, "../../figures/henon-heiles/periodic-orbits/")
const SAVE_DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/periodic-orbits/")

data_file        = joinpath(DATA_DIR, "E0.1127-T10000.0-py0.0-n256.jld2")


# ============================================================================================
#                       constants for precission of the simulations
# ============================================================================================


const CC_TOL  = 1e-13
const INT_TOL = 1e-14          # how precise the integrator should be for sattisfactory convergance
const PMAP_ROOT_TOL = 1e-11    # poincare returnmap precission to find periodic orbit, i.e. the root
const PMAP_PRIME_TOL = 1e-9    # tollerace for saying that two roots belong to the same orbit
const DEL_BOX = 1e-3           # offset between the two boxes that breaks the symmetry of the charged point particles


p = (;C=1, m1=1.0, m2=1.0, L=1.0, del= 1DEL_BOX)

E  = 1
v2 = [0.1,0.1]

u0 = lift(v2,E,p)
println("@$u0: E= ",energy(u0, p))


# params = SectionParams(E, p; 
#           tmax = 20_000.0, nfast = 1, ndense=40, 
#           cc_tol = CC_TOL, int_tol = INT_TOL, save_everystep = true, save_start = true)
int_time = 100_000
u, t, pts= get_traj(u0, int_time;
                    p=p, abstol=INT_TOL, reltol=INT_TOL)



xs = first.(u,2)
ps = last.(u,2)
xone= Point2f.(zip(t, first.(xs)))
xtwo= Point2f.(zip(t, last.(xs)))
sec = Point2f.(first.(pts,2))
sec_boundary  = boundary(E; p)
set_style!(:dark)
fig = Figure(size=(1300,900))

ax_sec = Axis(fig[2,1], xlabel=L"x_2\, [L]", ylabel=L"p_2")
ax_traj= Axis(fig[1,1:2],xlabel="time", ylabel="position in characteristic box size L", )
ax_energy= Axis(fig[2,2], xlabel="time", ylabel=L"\text{Energy}\,[\mathcal{C}/L]")
lines!(ax_traj, xone[1:2*2*500])
lines!(ax_traj, xtwo[1:2*2*500])
scatter!(ax_sec, sec)
scatter!(ax_sec, sec_boundary, color= INK)
lines!(ax_energy, t, [abs(energy(ui,p)- E)*p.L/E/p.C for ui in u] )
display(fig)