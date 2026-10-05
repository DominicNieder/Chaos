include("monodromy_method.jl")

energy_sweep = ABC_energy_trace()
orbits = append_monodrome(energy_sweep.all_ABC)

save_path = joinpath(@__DIR__, "../../data/henon-heiles/bifurcation/orbitsABC_E0.01-10_monodromy.jld2")
mkpath(dirname(save_path))
JLD2.save_object(save_path, orbits)
@info "Saved periodic orbits with monodromy matrices" path=save_path rows=nrow(orbits)


# ==========================
#       tori plot
# ==========================


set_style!(:print)  # :print
to = tori(orbits; labels=["A", "B", "C"],           n_periods=0, n_unst_perido=0, n_per_shell=8, n_shells=3, shell_radius = 3e-3, azimuth = pi/3, elevation = 0.9, perspectiveness = 0.3, mask=false)
display(to.fig)

# ==========================
#       orbits in config space
# ==========================
plt = show_orbits_In_Config(orbits; style=:print, n_trajectories=40)

display(plt.A.f)
display(plt.B.f)
display(plt.C.f)

# ==========================
#       following orbits and monodromy analysis
# ==========================
figs = graphs(orbits, save_as_svg = true)


display(figs.fC.fλ)
# display(figs.fB.fλ)
# display(figs.fC.fλ)

display(figs.fA.fLy)
# display(figs.fB.fLy)
# display(figs.fC.fLy)

display(figs.fC.fΘ)
# display(figs.fB.fΘ)
# display(figs.fC.fΘ)


# xlims!(figs.fB.ax, -0.001, 0.1)
# ylims!(figs.fB.ax, -0.5,3)
# save(joinpath(FIG_DIR, "orbitA/Aeigenvalues-vs-E.png"), figs.fA.fig; px_per_unit = 2)
