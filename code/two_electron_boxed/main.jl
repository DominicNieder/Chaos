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
const DEL_BOX = 1e-8          # offset between the two boxes that breaks the symmetry of the charged point particles





p =(;C=-1, m1=1.0, m2=1.0, L1=1.0, L2=1.0, del= 1e-8)

Es = collect(range(-0.5,-0.3,length= 256))
maps = scan_periodic_orbits_interactive(Es, p; nmax=3, nx=6, np=6,
                                          plims=(-5,5), background=true,
                                          save_fig=true, save_data=false,
                                          folder=joinpath(FIG_DIR, "periodic_orbits/scan2"),
                                          data_dir=DATA_DIR)

display(maps.fig)



# follow periodic orbits and analyse the monodromy matrix
Esfollow = sort!(collect(range(-2,-0.5, length=1000)),rev=true)
orbs = follow!(Esfollow, p, 3; seed= [0.113022,0.568])

ana = monodrome(orbs; p=p, verbose=true)
orbits_with_mono = copy(orbs)
orbits_with_mono.M = ana.M
orbits_with_mono.check = ana.check
orbits_with_mono.str .= "A"
set_style!(:dark)
eigenvalue_plots = graphs(orbits_with_mono; fsize=set_fig_size(0.8))
display(eigenvalue_plots.fA.fλ)

sec= section_map(-0.43, p; nx=5, np=5, tint=1_000)
display(sec.sec.f)





orb = maps.orbits
display(maps.maps[1].f)

seed = orb[(orb.E .== -1.0) .& ( orb.prime .== 1), :v]

EN1 = -0.42 /abs(p.C/p.L1)

prm = SectionParams(EN1, p; 
          tmax = 20_000.0, nfast = 1, ndense=40, 
          cc_tol = CC_TOL, int_tol = INT_TOL, save_everystep = false, save_start = false)
orb = analyse_seed(seedN1, 1, prm; str = "seed")
println(orb.v)
Es = collect(range(-20.0, 100.0, length=50000))

traj, _, _ = get_traj(lift(orb.v, EN1, p), orb.T; p=p)
delta = 1e-8
dtraj1, _, _ = get_traj(lift(orb.v, EN1, p) + [-delta, 0, 0, 0], orb.T; p=p)
dtraj2, _, _ = get_traj(lift(orb.v, EN1, p) + [0, delta, 0, 0], orb.T; p=p)
dtraj3, _, _ = get_traj(lift(orb.v, EN1, p) + [0, 0, delta, 0], orb.T; p=p)
dtraj4, _, _ = get_traj(lift(orb.v, EN1, p) + [0, 0, 0, delta], orb.T; p=p)

dtraj2p, _, _ = get_traj(lift(orb.v, EN1, p) + [0, -delta, 0, 0], orb.T; p=p)
dtraj3p, _, _ = get_traj(lift(orb.v, EN1, p) + [0, 0, -delta, 0], orb.T; p=p)
dtraj4p, _, _ = get_traj(lift(orb.v, EN1, p) + [0, 0, 0, -delta], orb.T; p=p)


uT = traj[end]
M = hcat(
    (dtraj1[end] .- uT) ./ (-delta),
    (dtraj2[end] .- dtraj2p[end]) ./ 2delta,
    (dtraj3[end] .- dtraj3p[end]) ./ 2delta,
    (dtraj4[end] .- dtraj4p[end]) ./ 2delta,
)
@show M det(M)
@show eigvals(M)

uen1 = flow(lift(orb.v, EN1, p) + [1e-8, 0, 0, 0], orb.T; p=p)

set_style!(:dark)
f = Figure(size=(900, 900))
ax = Axis(f[1, 1], xlabel="x1", ylabel="x2", aspect=DataAspect())
lines!(ax, Point2f.(getindex.(traj, 1), getindex.(traj, 2));
       color=INK, label="Unperturbed")
lines!(ax, Point2f.(getindex.(dtraj1, 1), getindex.(dtraj1, 2));
       color=pick_color(1), label="δx₁ = -1e-8")
scatter!(ax, [uen1[1]], [uen1[2]];
         markersize=8, color=pick_color(1), label="Endpoint (δx₁)")
lines!(ax, Point2f.(getindex.(dtraj2, 1), getindex.(dtraj2, 2));
       color=pick_color(2), label="δx₂ = 1e-8")
lines!(ax, Point2f.(getindex.(dtraj3, 1), getindex.(dtraj3, 2));
       color=pick_color(3), label="δp₁ = 1e-8")
lines!(ax, Point2f.(getindex.(dtraj4, 1), getindex.(dtraj4, 2));
       color=pick_color(4), label="δp₂ = 1e-8")
axislegend(ax; position=:rt)

display(f)


# so I want to find the periodic orbit 

# ================================
# for scanning the section maps
# ================================
# scan = scan_section_maps(Es ./ abs(p.C/p.L1), p; nx=6, np=2, tint=5_000, save_fig=true, plims=(-40,40))

# display(scan[1].sec.f)
# display(scan[2].sec.f)
# display(scan[3].sec.f)
# display(scan[4].sec.f)
# display(scan[5].sec.f)
# display(scan[6].sec.f)
# display(scan[7].sec.f)
# display(scan[8].sec.f)
# display(scan[9].sec.f)
# display(scan[10].sec.f)
# diplay(sca)