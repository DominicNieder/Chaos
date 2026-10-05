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

del = [0.0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5]
p = [(;C=-1, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= deli) for deli in del]


int_time = 10_000

res = explore_del_between_boxes(p;style=:print, save_fig=true, int_time=1000)


# display(res[1].f)