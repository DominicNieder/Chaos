# =====================================================================
#  Driver for the periodic-orbit sweep.
#
#  `symmetries.jl` holds definitions only and has no side effects beyond the
#  load-time asserts, so it can be included by tests, by notebooks, or here.
#  This file is the only place a sweep is actually started.
#
#      julia --project=.. code/finding_period_orbits/run_periodic_orbits.jl
#
#  or, interactively:
#
#      include("run_periodic_orbits.jl")
# =====================================================================

include("test.jl")
include("a_tester.jl")

# @testset throws a TestSetException on failure, so a broken invariant stops
# the run here rather than after several hours of integration.
run_tests()

const RUN_SLOW_TESTS = false      # a few seconds; integrates the flow
RUN_SLOW_TESTS && run_tests_slow(0.12)


# ---------------------------------------------------------------------
#  The sweep
# ---------------------------------------------------------------------

# res = main(
#     Emax        = 0.166666,
#     Emin        = 5e-2,
#     n_energies  = 166 * 10,
#     ns          = 1:9,
#     tmax        = 20_000.0,
#     match       = :exact,
#     reduce_to   = :classes,
#     outfile     = "05Sept_test_following_orbits_classes.jld2",
#     save_data   = true,
# )


# ---------------------------------------------------------------------
#  Quick look at the result
# ---------------------------------------------------------------------
f = joinpath(SAVE_DATA_DIR, "05Sept_test_following_orbits_classes.jld2")   # orbits_E0.16-0.0001_160.jld2
orbits, timing, energies, reduce_to, match, group_order, newton_tol, ode_abstol = load(f, "orbits", "timing", "energies", "reduce_to", "match", "group_order", "newton_tol", "ode_abstol")

const E_TOP = maximum(energies)

traced = trace_branches(orbits)
println(first(branch_report(traced), 20))
tra
TE = plot_T_vs_E(traced; min_len = 3)
display(TE.fig)

prm  = SectionParams(E_TOP, PARAM; tmax = 20_000.0, ndense = NMAX_SEARCH)
top  = at_energy(orbits, E_TOP)
full = expand_symmetry(top, prm)

display(plot_orbits(full, prm).fig)
display(section_slider(orbits, PARAM).fig)