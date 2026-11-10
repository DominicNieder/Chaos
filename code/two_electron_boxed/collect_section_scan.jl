include("two_electron_boxed_apply_monodromy.jl")
include("analysis.jl")

BLAS.set_num_threads(1)          # your linear algebra is 2x2; BLAS threads only compete


# ============================================================================================
#                       Save and data directories
# ============================================================================================

const DATA_DIR   = joinpath(@__DIR__, "../../data/two-boxed-charges/pSec/")

const INTERVALL1 = joinpath(DATA_DIR, "-infinityto-1/")
const INTERVALL2 = joinpath(DATA_DIR, "-1to0/")
const INTERVALL3 = joinpath(DATA_DIR, "0to+infinity/")

mkpath(INTERVALL1)
mkpath(INTERVALL2)
mkpath(INTERVALL3)

# ============================================================================================
#                       constants for precission of the simulations
# ============================================================================================


const CC_TOL  = 1e-13
const INT_TOL = 1e-14          # how precise the integrator should be for sattisfactory convergance
const PMAP_ROOT_TOL = 1e-11    # poincare returnmap precission to find periodic orbit, i.e. the root
const PMAP_PRIME_TOL = 1e-9    # tollerace for saying that two roots belong to the same orbit
const DEL_BOX = 1e-8          # offset between the two boxes that breaks the symmetry of the charged point particles



nx, np = 6, 10 

p =(;C=-1, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= DEL_BOX)

number_of_slices = 128

n1 = number_of_slices
n2 = number_of_slices
n3 = number_of_slices
n4 = number_of_slices

tinit = 2000

# E_scan1 = collect(range(-500, -1.0, length=n1))
# E_scan2 = collect(range(-2, -0.3, length=n2))
E_scan3 = collect(range(-0.3, 2, length=n3))
E_scan4 = collect(range(2, 5_000, length=n4))

scan_dir = joinpath(FIG_DIR, "section_maps/scanAll2")
mkpath(scan_dir)

energy_ranges = (E_scan3, E_scan4)  
E_scan = vcat(energy_ranges...)

# Simulate all four ranges in one sequence so the data and video share one slider.
sec_scan = interactive_section_maps(E_scan, p; nx=nx, np=np, tint=tinit, seeds=:relative,
                                    plims=nothing, c=0.999, vis=:dark,
                                    markersize=3, parallel=false, show_progress=true,
                                    folder=scan_dir, drift_on=false)


display(sec_scan.fig)

scan_data = joinpath(scan_dir, "section_scan.jld2")
jldsave(scan_data; energies=sec_scan.energies, energy_ranges, parameters=p,
        seeds=sec_scan.seeds, sections=sec_scan.secs, drift=sec_scan.drift,
        settings=(; nx, np, tint=tinit, seeds=:relative, c=0.999,
                  plims=nothing, vis=:dark, markersize=3, drift_on=false))
@info "Saved section-map scan data" scan_data

scan_video = joinpath(scan_dir, "section_scan.mp4")
record(sec_scan.fig, scan_video, eachindex(sec_scan.energies); framerate=10) do i
    set_close_to!(sec_scan.slider, i)
end
@info "Saved section-map scan video" scan_video
