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

# Es = collect(range(-1.5,-0.3,length= 256))
# secs = interactive_section_maps(Es, p; 
#                                    nx=20, np=20,              tint=1_000,           seeds=:relative,
#                                   plims=(-5,5), follow_limits=true, vis=:dark,
#                                   markersize=4, parallel=true, show_progress=true,
#                                   folder=joinpath(FIG_DIR, "section_maps/scan1/section_scan.mp4"))
# display(secs.fig)
# mkpath(joinpath(FIG_DIR, "section_maps/scan1"))
# record(secs.fig, joinpath(FIG_DIR, "section_maps/scan1/section_scan.mp4"), eachindex(secs.energies); framerate = 10) do i
#     set_close_to!(secs.slider, i)
# end

# ========================================
#  following specific periodic orbits
# ========================================
# defining 'following' parameters
EsA1 = collect(range(-1.1,-0.3, length=1000))
labelA1, seedA1, periodA1 = "A1", [0.02,0.0], 1

EsA31 = sort!(collect(range(-0.5, -0.0, length=1000)), rev=true)
labelA3, seedA3, periodA3 = "A3.1", [0.118, 0.57], 3

# following specified parameters
# A1
orbsA1 = follow!(EsA1, p, periodA1; 
       seed= seedA1, 
       str= labelA1)
# A3.1
orbsA31 = follow!(EsA31, p, periodA3; 
       seed= seedA3, 
       str= labelA3)





# pmap= section_map(e, p; nx=20, np=20, tint=1_000)

# mkpath(joinpath(FIG_DIR, "section_maps/A3.1"))
# save(joinpath(FIG_DIR, "section_maps/A3.1/E-$(round(e,5)).png"), pmap.sec.f)

# emax  = maximum(orbitA31_plot.orbits.E[orbitA31_plot.orbits.jump]) + 0.001
# emin  = minimum(orbitA31_plot.orbits.E[orbitA31_plot.orbits.jump])- 0.001
# es = collect(range(emin, emax, 64))
# int_pmap = interactive_section_maps(es, p; nx=20, np=20, tint=1_000, seeds=:relative,
#                                   plims=nothing, follow_limits=true, vis=:dark,
#                                   markersize=3, parallel=true, show_progress=true,
#                                   folder=joinpath(FIG_DIR, "section_maps/A3.1"))

# record(int_pmap.fig, joinpath(FIG_DIR, "section_maps/A3.1/creation_zoom.mp4"), eachindex(int_pmap.energies); framerate = 10) do i
#     set_close_to!(int_pmap.slider, i)
#     limits!(int_pmap.ax, 0.0, 0.2, 0, 1)          # (xmin, xmax, ymin, ymax), fixed for all frames
# end
# on(int_pmap.slider.value) do i
#     limits!(int_pmap.ax, 0.0, 1, -5, 5)              # fixed frame
# end
# notify(int_pmap.slider.value)                          # apply now, not only on the next step
# display(int_pmap.fig)

orb = maps.orbits
display(maps.maps[1].f)



orbitA1_plot = analyse_followed_traj(orbsA1, p; trajectory_every=50)
display(orbitA1_plot.f)
display(orbitA1_plot.config_f)
display(getproperty(orbitA1_plot.monodrome, Symbol("f", first(orbsA1.str))).fΘ)

orbitA31_plot = analyse_followed_traj(orbsA31, p; trajectory_every=50)
display(orbitA31_plot.f)
display(orbitA31_plot.config_f)
display(getproperty(orbitA31_plot.monodrome, Symbol("f", first(orbsA31.str))).fλ)


# studying the orbit at creation

orbits = orbitA31_plot.orbits

jump_indices = findall(orbits.jump)
jump_orb= orbits[jump_indices,:]
idx = jump_indices[argmin(orbits.E[jump_indices])]
println(jump_orb)
orbit = orbits[idx, :]
println("Index: ", idx)
println("Energy: ", orbit.E)
println("Initial condition: ", orbit.v)
println("Initial theta: ", orbit.theta)
emin  = minimum(orbitA31_plot.orbits.E[orbitA31_plot.orbits.jump])
emax  = maximum(orbitA31_plot.orbits.E[orbitA31_plot.orbits.jump])
println(emin, emax)
es = collect(range(emin, emax, 128))

s = section_map(jump_orb.E[end-6], p; nx=10, np=10, tint=250, plims=(-5,5), set_xmin=0.04, set_xmax=0.06, set_pmin=0.45, set_pmax= 0.48)
display(s.sec.f)
mkpath(joinpath(FIG_DIR, "section_maps_attractive/detailE-051"))
save(joinpath(FIG_DIR, "section_maps_attractive/detailE-051/sec.png"),s.sec.f)



A31_detail = follow!(es,p,3, seed=orbit.v)


num = 3
seed = orbsA31.v[num]
EN1 = orbsA31.E[num]/abs(p.C/p.L1)

u, t, pts = get_traj(lift(seed,EN1,p), orbsA31.T[num];p=p)
set_style!(:dark)
f = Figure(size=(900, 900))
ax = Axis(f[1, 1], xlabel="x1", ylabel="x2", aspect=DataAspect())

lines!(ax, Point2f.(getindex.(u, 1).* sign.(getindex.(u,4)), getindex.(u, 2)) .* sign.(getindex.(u,4));
       color=INK, label="Unperturbed")
display(f)

prm = SectionParams(EN1, p; 
          tmax = 20_000.0, nfast = 1, ndense=40, 
          cc_tol = CC_TOL, int_tol = INT_TOL, save_everystep = false, save_start = false)
orb = analyse_seed(seed, 1, prm; str = "seed")
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