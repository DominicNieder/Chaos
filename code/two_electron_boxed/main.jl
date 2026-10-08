include("two_electron_boxed_apply_monodromy.jl")
include("analysis.jl")

BLAS.set_num_threads(1)          # your linear algebra is 2x2; BLAS threads only compete


# ============================================================================================
#                       Save and data directories
# ============================================================================================

const DATA_DIR   = joinpath(@__DIR__, "../../data/two-boxed-charges/")



# ============================================================================================
#                       constants for precission of the simulations
# ============================================================================================


const CC_TOL  = 1e-13
const INT_TOL = 1e-14          # how precise the integrator should be for sattisfactory convergance
const PMAP_ROOT_TOL = 1e-11    # poincare returnmap precission to find periodic orbit, i.e. the root
const PMAP_PRIME_TOL = 1e-9    # tollerace for saying that two roots belong to the same orbit
const DEL_BOX = 1e-15          # offset between the two boxes that breaks the symmetry of the charged point particles

p =(;C=-1, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= 1e-8)

Es = [-0.55,-0.45,-0.4, -0.3, -0.2] # [-500, -100, -50, -10, -5,-1,-0.1,0.1,0.5, 1.0, 2.0, 3.0, 4.0]
set_style!(:dark)



# po = scan_periodic_orbits(Es ./ abs(p.C/p.L1), p; nmax=8, nx=5, np=5, save_fig=true,plims=(-20,20))
# display(po.maps[4].f)



scan = scan_section_maps(Es ./ abs(p.C/p.L1), p; nx=6, np=2, tint=5_000, save_fig=true, plims=(-40,40))

display(scan[1].sec.f)
display(scan[2].sec.f)
display(scan[3].sec.f)
display(scan[4].sec.f)
display(scan[5].sec.f)
display(scan[6].sec.f)
display(scan[7].sec.f)
display(scan[8].sec.f)
display(scan[9].sec.f)
display(scan[10].sec.f)
# diplay(sca)