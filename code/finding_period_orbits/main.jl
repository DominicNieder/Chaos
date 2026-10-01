include("monodromy_method.jl")


# res = ABC_energy_trace(nup=5000, ndown=5000).all_ABC
# orbits = append_monodrome(res)



set_style!(:print)  # :print
to = tori(orbits; labels=["A", "B", "C"],           n_periods=0, n_unst_perido=0, n_per_shell=8, n_shells=3, shell_radius = 3e-3, azimuth = pi/3, elevation = 0.9, perspectiveness = 0.3, mask=false)
display(to.fig)



figs = graphs(orbits)


display(figs.fA.fλ)
# display(figs.fB.fλ)
# display(figs.fC.fλ)

display(figs.fA.fLy)
# display(figs.fB.fLy)
# display(figs.fC.fLy)

display(figs.fA.fΘ)
# display(figs.fB.fΘ)
# display(figs.fC.fΘ)


# xlims!(figs.fB.ax, -0.001, 0.1)
# ylims!(figs.fB.ax, -0.5,3)
# save(joinpath(FIG_DIR, "orbitA/Aeigenvalues-vs-E.png"), figs.fA.fig; px_per_unit = 2)
