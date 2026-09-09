# =====================================================================
#  Tests for symmetries.jl
#
#  Include this AFTER the definitions and BEFORE the driver:
#
#      include("symmetries.jl")     # with the `res = main(...)` line commented out
#      include("test_symmetries.jl")
#      run_tests()
#
#  Everything in `run_tests()` is pure logic: no ODE integration, no Newton
#  solves, runs in well under a second. The integration-level checks that do
#  need the flow are in `run_tests_slow(E)`, kept separate so the fast set can
#  be run on every edit.
# =====================================================================

using Test, LinearAlgebra, DataFrames, Random


# ---------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------

"Total energy of a 4D state, for invariance checks."
test_H(u, p = PARAM) = (u[3]^2 + u[4]^2) / (2p[2]) +
                       HenonHeiles.potential(u[1], u[2], p)

"4x4 matrix of the phase-space action of g, in (q, p) ordering."
function test_action_matrix(g)
    A = zeros(4, 4)
    for k in 1:4
        e = zeros(4); e[k] = 1.0
        A[:, k] = apply_sym(g, e)
    end
    return A
end

"The symplectic form in (q1, q2, p1, p2) ordering. A function, not a const, so
this file can be re-included in a live session without redefinition warnings."
test_J() = [zeros(2, 2) I(2); -I(2) zeros(2, 2)]

"""
A synthetic `Row`. Pushing one into `orbit_table()` is itself a schema test:
it fails if `Row`'s field names, order or types drift away from what
`analyse_seed` produces.
"""
function test_row(; y, py, sec_y = [y], sec_py = [py],
                    sym_y = Float64[], sym_py = Float64[],
                    prime = 1, T = 2.0, trace = 0.5, E = 0.1,
                    resnorm = 1e-10, index = 1, sym_stab = "", sym_mult = 0,
                    sym_id = 0, id = 0, origin = "test")
    return (; E = E, n = prime, seed_y = y, seed_py = py,
              y = y, py = py, prime = prime, T = T,
              trace = trace, detDT = 1.0,
              resnorm = resnorm, iters = 1,
              sec_y = sec_y, sec_py = sec_py,
              history_y = [y], history_py = [py],
              traj_x = [0.0], traj_y = [y],
              class = KIND_LABEL[index], index = index,
              sym_y = sym_y, sym_py = sym_py,
              sym_stab = sym_stab, sym_mult = sym_mult, sym_id = sym_id,
              id = id, origin = origin)
end

test_frame(rows...) = DataFrame(collect(rows))


# ---------------------------------------------------------------------
#  Fast tests
# ---------------------------------------------------------------------

function run_tests()
@testset "symmetries.jl" begin

# =====================================================================
@testset "group structure" begin
    @test length(SYM_GROUP) == 12
    @test length(unique(g.name for g in SYM_GROUP)) == 12

    # every element is an orthogonal configuration-space map
    for g in SYM_GROUP
        @test norm(g.M' * g.M - I) < 1e-12
        @test abs(abs(det(g.M)) - 1) < 1e-12
        @test g.τ in (1, -1)
    end

    # closure, and τ is multiplicative (the Z2 factor is a direct product)
    for g in SYM_GROUP, h in SYM_GROUP
        k = compose(g, h)
        @test k in SYM_GROUP
        @test k.τ == g.τ * h.τ
        @test norm(k.M - g.M * h.M) < 1e-9
    end

    # identity and inverses
    e = only(filter(g -> norm(g.M - I) < 1e-12 && g.τ == 1, SYM_GROUP))
    @test e.name == "E"
    for g in SYM_GROUP
        @test any(h -> compose(g, h).name == "E", SYM_GROUP)
    end

    # Lagrange: the order of every cyclic subgroup divides |G| = 12
    for g in SYM_GROUP
        k, ord = g, 1
        while k.name != "E" && ord <= 12
            k = compose(k, g); ord += 1
        end
        @test ord <= 12
        @test 12 % ord == 0
    end

    # D3 has three rotations and three reflections
    rotations = count(g -> g.τ == 1 && det(g.M) > 0, SYM_GROUP)
    @test rotations == 3
    @test count(g -> g.τ == 1 && det(g.M) < 0, SYM_GROUP) == 3

    # C3 really is order 3
    c3 = only(filter(g -> g.name == "C3", SYM_GROUP))
    @test compose(compose(c3, c3), c3).name == "E"
end

# =====================================================================
@testset "group action preserves the dynamics" begin
    rng = MersenneTwister(20260905)

    for g in SYM_GROUP, _ in 1:50
        u = 0.6 .* (2 .* rand(rng, 4) .- 1)
        @test isapprox(test_H(apply_sym(g, u)), test_H(u); atol = 1e-12)
    end

    # A^T J A = τ J: symplectic when τ = +1, anti-symplectic when τ = -1.
    J = test_J()
    for g in SYM_GROUP
        A = test_action_matrix(g)
        @test norm(A' * J * A - g.τ * J) < 1e-12
    end

    # the potential is invariant (this is what check_group asserts at load)
    @test check_group() === true

    # y -> -y alone is NOT a symmetry: the guard against the obvious mistake
    bad = [1.0 0.0; 0.0 -1.0]
    @test !isapprox(HenonHeiles.potential((bad * [0.3, 0.4])..., PARAM),
                    HenonHeiles.potential(0.3, 0.4, PARAM); atol = 1e-8)
end

# =====================================================================
@testset "section stabiliser and cosets" begin
    @test length(SEC_STAB) == 2
    @test Set(g.name for g in SEC_STAB) == Set(["E", "σΘ"])
    @test length(COSET_REPS) == 6

    # the 6 cosets tile the group exactly: 6 x 2 = 12 distinct elements
    tiled = [compose(g, h).name for g in COSET_REPS for h in SEC_STAB]
    @test length(tiled) == 12
    @test length(unique(tiled)) == 12
    @test Set(tiled) == Set(g.name for g in SYM_GROUP)

    # only E and σΘ pass the predicate
    @test count(is_section_stabiliser, SYM_GROUP) == 2

    # σΘ acts as (y, py) -> (y, -py) and preserves px^2, so the partner point
    # is automatically inside the energy boundary
    s = only(filter(g -> g.name == "σΘ", SEC_STAB))
    v = [0.3, 0.2]
    @test sec_action(s, v) ≈ [0.3, -0.2]
    E = 0.12
    @test px2(v[1], v[2], E, PARAM) ≈ px2(sec_action(s, v)..., E, PARAM)
    @test in_section(v, E, PARAM) == in_section(sec_action(s, v), E, PARAM)
end

# =====================================================================
@testset "section geometry" begin
    E, p = 0.12, PARAM

    # px^2 is what is left of the energy at x = 0
    y, py = 0.2, 0.1
    @test px2(y, py, E, p) ≈ 2p[2] * (E - HenonHeiles.potential(0.0, y, p)) - py^2

    # the lift lands on the energy shell (up to the O(eps^2) x-offset)
    v = [0.2, 0.1]
    u = lift(v, E, p)
    @test u !== nothing
    @test isapprox(test_H(u, p), E; atol = 1e-12)
    @test u[1] > 0 && u[3] > 0          # x and px both positive: leaves the section
    @test u[2] == v[1] && u[4] == v[2]

    # sgn = -1 flips both, so the state again moves away from x = 0
    um = lift(v, E, p; sgn = -1)
    @test um[1] < 0 && um[3] < 0
    @test isapprox(test_H(um, p), E; atol = 1e-12)

    # outside the boundary
    @test lift([0.2, 10.0], E, p) === nothing
    @test !in_section([0.2, 10.0], E, p)

    # pymax is exactly where px^2 vanishes
    pm = pymax(0.2, E, p)
    @test isapprox(px2(0.2, pm, E, p), 0.0; atol = 1e-12)
    @test in_section([0.2, 0.99pm], E, p)
    @test !in_section([0.2, 1.01pm], E, p)

    # margin tightens the test
    @test !in_section([0.2, 0.999999pm], E, p; margin = 1e-4)
end

# =====================================================================
@testset "stability classification" begin
    @test kind_index(0.0)  == 1     # elliptic
    @test kind_index(1.5)  == 1
    @test kind_index(-1.5) == 1
    @test kind_index(3.0)  == 2     # hyperbolic
    @test kind_index(-3.0) == 3     # inverse hyperbolic
    @test kind_index(2.0)  == 4     # parabolic
    @test kind_index(-2.0) == 4
    @test kind_index(2.0 + 1e-9) == 4
    @test kind_index(2.1) == 2      # outside the parabolic window
    @test all(1 .<= kind_index.([0.0, 3.0, -3.0, 2.0]) .<= length(KIND_LABEL))
end

# =====================================================================
@testset "symmetry type names" begin
    @test sym_type("") == "unknown"
    @test sym_type(join([g.name for g in SYM_GROUP], ",")) == "fully symmetric"
    @test sym_type("E,C3,C3²") == "C3 invariant"
    @test sym_type("E,σ") == "mirror invariant"
    @test sym_type("E,σΘ") == "reversible"
    @test sym_type("E,σ,Θ,σΘ") == "mirror + reversible"
    @test sym_type("E") == "generic"
end

# =====================================================================
@testset "Row schema" begin
    df = orbit_table()
    @test nrow(df) == 0
    # push! only succeeds if names, order and types all match Row
    @test_nowarn push!(df, test_row(y = 0.1, py = 0.2))
    @test nrow(df) == 1
    for f in (:E, :prime, :T, :trace, :sym_y, :sym_py, :sym_stab,
              :sym_mult, :sym_id, :id, :origin)
        @test f in propertynames(df)
    end
    @test eltype(df.prime) == Int
    @test eltype(df.sym_y) == Vector{Float64}
end

# =====================================================================
@testset "deduplication" begin
    a = test_row(y = 0.1, py = 0.2)
    b = test_row(y = 0.1 + 1e-8, py = 0.2 - 1e-8)

    @testset "literal duplicates" begin
        @test same_point_set(a, b; tol = 1e-6)
        @test !same_point_set(a, test_row(y = 0.5, py = 0.2); tol = 1e-6)
        @test nrow(dedup_raw(test_frame(a, b); tol = 1e-6)) == 1
        # different prime periods are never merged, however close the points
        @test nrow(dedup_raw(test_frame(a, test_row(y = 0.1, py = 0.2, prime = 2));
                             tol = 1e-6)) == 2
    end

    @testset "symmetry partners" begin
        pa = test_row(y = 0.3, py = 0.4, sym_y = [0.7], sym_py = [-0.2])
        pb = test_row(y = 0.7, py = -0.2, T = 2.01, trace = 0.7)
        @test symmetry_match(pa, pb; tol = 1e-6)
        @test symmetry_match(pb, pa; tol = 1e-6)      # must be symmetric
        @test same_class(pa, pb; tol = 1e-6)
        @test same_class(pb, pa; tol = 1e-6)
        # annotation on only one side is enough
        @test same_class(test_row(y = 0.7, py = -0.2), pa; tol = 1e-6)

        @test !same_class(pa, test_row(y = 0.8, py = -0.1); tol = 1e-6)
        @test !same_class(pa, test_row(y = 0.7, py = -0.2, prime = 2); tol = 1e-6)
        @test !same_class(pa, test_row(y = 0.7, py = -0.2, E = 0.2); tol = 1e-6)

        d = class_diagnostics(pa, pb)
        @test d.same_geometry && d.same_prime && d.same_energy
        @test !d.same_period && !d.same_trace
    end

    @testset "order independence" begin
        rng = MersenneTwister(7)
        rows = [test_row(y = 0.1, py = 0.0),
                test_row(y = 0.1 + 1e-9, py = 0.0),
                test_row(y = 0.5, py = 0.0),
                test_row(y = 0.5 + 1e-9, py = 0.0),
                test_row(y = 0.9, py = 0.0)]
        counts = [nrow(dedup_raw(DataFrame(shuffle(rng, rows)); tol = 1e-6))
                  for _ in 1:20]
        @test all(==(3), counts)      # union-find must not depend on row order
    end

    @testset "class representatives" begin
        pa = test_row(y = 0.3, py = 0.4, sym_y = [0.7], sym_py = [-0.2],
                      resnorm = 1e-9)
        pb = test_row(y = 0.7, py = -0.2, resnorm = 1e-11)
        pc = test_row(y = -0.4, py = 0.2, resnorm = 1e-8)
        cls = dedup_classes(test_frame(pa, pb, pc); tol = 1e-6)
        @test nrow(cls) == 2
        @test cls.sym_id == [1, 2]
        # sorted by resnorm, so the best-converged member represents the class
        @test cls.resnorm[1] ≈ 1e-11
    end

    @testset "energies are never merged" begin
        lo = test_row(y = 0.3, py = 0.4, E = 0.10)
        hi = test_row(y = 0.3, py = 0.4, E = 0.11)
        @test !same_class(lo, hi; tol = 1e-6)
        @test nrow(dedup_all(test_frame(lo, hi))) == 2
    end

    @testset "transitive closure is a tolerance risk" begin
        # Union-find merges A~B and B~C into one component even when A and C are
        # further apart than tol. This is correct for an exact equivalence, but
        # with an approximate one it can chain across a near-collision of two
        # genuinely distinct orbits (as happens near a tangent bifurcation).
        # The test pins the current behaviour so a design change is visible.
        chain = test_frame(test_row(y = 0.0, py = 0.0),
                           test_row(y = 0.1, py = 0.0),
                           test_row(y = 0.2, py = 0.0))
        @test nrow(dedup_by(chain, (x, z) -> abs(x.y - z.y) <= 0.11)) == 1
        @test nrow(dedup_by(chain, (x, z) -> abs(x.y - z.y) <= 0.05)) == 3
    end

    @testset "no spurious merging" begin
        rows = [test_row(y = 0.1i, py = 0.0) for i in 1:6]
        @test nrow(dedup_raw(DataFrame(rows); tol = 1e-6)) == 6
        @test nrow(dedup_classes(DataFrame(rows); tol = 1e-6)) == 6
    end
end

# =====================================================================
@testset "energy grids" begin
    Es = log_energies(1e-4, 1e-1; per_decade = 10)
    @test length(Es) == 31                      # 3 decades x 10 + 1
    @test issorted(Es)
    @test Es[1] ≈ 1e-4 && Es[end] ≈ 1e-1

    ρ = Es[2:end] ./ Es[1:end-1]
    @test maximum(ρ) - minimum(ρ) < 1e-9        # constant ratio ...
    d = diff(log10.(Es))
    @test maximum(d) - minimum(d) < 1e-12       # ... i.e. constant Δlog E

    @test length(log_energies(1e-4, 1e-1; n = 7)) == 7
    @test_throws ArgumentError log_energies(1e-1, 1e-4; n = 5)      # reversed
    @test_throws ArgumentError log_energies(0.0, 1e-1; n = 5)       # Emin = 0
    @test_throws ArgumentError log_energies(1e-4, 1e-1)             # neither
    @test_throws ArgumentError log_energies(1e-4, 1e-1; n = 5, per_decade = 10)
    @test_throws ArgumentError log_energies(0.16664, 0.16665; per_decade = 200)

    Ec = 1/6
    Ecs = crit_energies(Ec, 1e-6, 1e-3; per_decade = 10)
    @test issorted(Ecs)
    @test all(Ecs .< Ec)
    dd = diff(log10.(Ec .- Ecs))
    @test maximum(dd) - minimum(dd) < 1e-9      # geometric in the distance to Ec

    r = grid_report(Es)
    @test r.ratio[2] - r.ratio[1] < 1e-9
end

# =====================================================================
@testset "seed routing" begin
    pts    = [[0.1, 0.0], [0.2, 0.0], [0.3, 0.0], [0.4, 0.0]]
    primes = [1, 2, 2, 4]
    tagged = collect(zip(pts, primes))

    ex = seed_router(tagged, nothing, 0.1, PARAM; match = :exact)
    @test length(ex(1)) == 1
    @test length(ex(2)) == 2
    @test length(ex(3)) == 0
    @test length(ex(4)) == 1

    dv = seed_router(tagged, nothing, 0.1, PARAM; match = :divisors)
    @test length(dv(1)) == 1                    # 1 | 1
    @test length(dv(2)) == 3                    # 1 | 2, 2 | 2
    @test length(dv(4)) == 4                    # 1, 2, 4 all divide 4
    @test length(dv(3)) == 1                    # only 1 | 3

    al = seed_router(tagged, nothing, 0.1, PARAM; match = :all)
    @test length(al(1)) == length(al(7)) == 4

    # parallel-vector form must route identically to the tagged form
    par = seed_router(pts, primes, 0.1, PARAM; match = :exact)
    @test length(par(2)) == 2
    @test par(2) == ex(2)

    # an untagged flat list is used for every n
    flat = seed_router(pts, nothing, 0.1, PARAM)
    @test length(flat(1)) == length(flat(5)) == 4

    @test_throws ArgumentError seed_router(pts, [1, 2], 0.1, PARAM)
    @test_throws ArgumentError seed_router(tagged, nothing, 0.1, PARAM;
                                           match = :nonsense)

    # the router never mutates its input
    @test primes == [1, 2, 2, 4]
end

# =====================================================================
@testset "row selection" begin
    df = test_frame(test_row(y = 0.1, py = 0.0, E = 0.10),
                    test_row(y = 0.2, py = 0.0, E = 0.10),
                    test_row(y = 0.3, py = 0.0, E = 0.11))
    @test nrow(at_energy(df, 0.10)) == 2
    @test nrow(at_energy(df, 0.11)) == 1
    @test_throws ErrorException at_energy(df, 0.5)
    # rtol, not exact equality: a float rebuilt from the grid still matches
    @test nrow(at_energy(df, 0.10 * (1 + 1e-12))) == 2
end

# =====================================================================
@testset "branch tracing" begin
    # one orbit drifting slowly with energy: a single branch
    d = trace_branches(test_frame(
            test_row(y = 0.30, py = 0.10, E = 0.11, T = 2.00),
            test_row(y = 0.31, py = 0.10, E = 0.10, T = 2.02),
            test_row(y = 0.32, py = 0.10, E = 0.09, T = 2.04)))
    @test length(unique(d.branch)) == 1
    @test "branch" in names(d)

    # a jump larger than tol starts a new branch
    d2 = trace_branches(test_frame(
            test_row(y = 0.30, py = 0.10, E = 0.11),
            test_row(y = 0.90, py = 0.10, E = 0.10)); tol = 0.05)
    @test length(unique(d2.branch)) == 2

    # different prime periods never link, however close the points
    d3 = trace_branches(test_frame(
            test_row(y = 0.30, py = 0.10, E = 0.11, prime = 1),
            test_row(y = 0.30, py = 0.10, E = 0.10, prime = 2)))
    @test length(unique(d3.branch)) == 2

    # symmetry-aware: dedup_classes picked a different class member at the
    # lower energy, but the predecessor's stored roots still match it
    d4 = trace_branches(test_frame(
            test_row(y = 0.30, py = 0.10, E = 0.11,
                     sym_y = [0.30, -0.55], sym_py = [0.10, 0.22]),
            test_row(y = -0.55, py = 0.22, E = 0.10)))
    @test length(unique(d4.branch)) == 1

    # two coexisting orbits stay two branches across three energies
    rows = vcat([[test_row(y = 0.30 + 0.01k, py = 0.10, E = 0.11 - 0.01k),
                  test_row(y = -0.50 - 0.01k, py = 0.20, E = 0.11 - 0.01k)]
                 for k in 0:2]...)
    d5 = trace_branches(DataFrame(rows))
    @test length(unique(d5.branch)) == 2
    br = branch_report(d5)
    @test all(br.nE .== 3)
end

end # @testset "symmetries.jl"
end # run_tests


# ---------------------------------------------------------------------
#  Slow tests: these integrate the flow
# ---------------------------------------------------------------------

"""
    run_tests_slow(E = 0.12; tmax = 5000.0)

Checks that need the ODE. A few seconds rather than milliseconds, so they are
not part of `run_tests()`.
"""
function run_tests_slow(E = 0.12; tmax = 5000.0)
@testset "integration-level" begin
    prm = SectionParams(E, PARAM; tmax, nfast = 1, ndense = 20)

    @testset "the return map returns to the section" begin
        v = [0.2, 0.05]
        w = poincare_map(v, 1, prm)
        @test length(w) == 2
        @test in_section(w, E, PARAM)
        # the image is on the same energy shell
        @test isapprox(test_H(lift(w, E, PARAM), PARAM), E; atol = 1e-10)
        # T^2 = T o T
        @test isapprox(poincare_map(v, 2, prm),
                       poincare_map(w, 1, prm); atol = 1e-7)
    end

    @testset "with_ndense restores state" begin
        before = prm.nmax_dense[]
        with_ndense(prm, 99) do
            @test prm.nmax_dense[] == 99
        end
        @test prm.nmax_dense[] == before
        # and restores it even when the body throws
        @test_throws ErrorException with_ndense(prm, 99) do
            error("boom")
        end
        @test prm.nmax_dense[] == before
    end

    @testset "a converged root really is periodic" begin
        seeds = section_grid(E, PARAM; ny = 4, npy = 4)
        got = nothing
        for v0 in seeds
            r = try solve_orbit(v0, 1, prm) catch; nothing end
            r === nothing && continue
            r.converged && (got = r; break)
        end
        if got === nothing
            @warn "no period-1 orbit converged at E = $E; skipping"
        else
            @test norm(Fres(got.v, 1, prm)) < 1e-8
            @test isapprox(det(got.DT), 1.0; atol = 1e-4)   # area preserving
        end
    end

    @testset "symmetry images are genuine orbits" begin
        seeds = section_grid(E, PARAM; ny = 4, npy = 4)
        row = nothing
        for v0 in seeds
            r = try analyse_seed(v0, 1, prm) catch; nothing end
            r === nothing && continue
            row = r; break
        end
        if row === nothing
            @warn "no orbit found at E = $E; skipping"
        else
            s = symmetry_data([row.y, row.py], row.sec_y, row.sec_py, prm)
            @test !isempty(s.stab)
            @test "E" in s.stab                       # identity always stabilises
            @test 12 % length(s.stab) == 0            # Lagrange
            @test s.mult in (1, 2, 3, 4, 6, 12)
            @test length(s.roots) <= 12
            # every image root is on the energy shell and inside the boundary
            for w in s.roots
                @test in_section(w, E, PARAM)
            end
            # each image is a period-n orbit with the same T and trace
            for w in s.roots
                img = try analyse_seed(w, row.prime, prm) catch; nothing end
                img === nothing && continue
                @test img.prime == row.prime
                @test isapprox(img.T, row.T; rtol = 1e-6)
                @test isapprox(img.trace, row.trace; rtol = 1e-3, atol = 1e-6)
            end
        end
    end
end
end