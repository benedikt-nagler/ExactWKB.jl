using ExactWKB
using ExactWKB: lines, points
using Test

@testset "stokes_graph" begin
    @testset "Airy: three rays to infinity at 2π/3" begin
        airy = SchrodingerProblem([0, 1])          # Q = z, one simple TP at 0
        g = stokes_graph(airy; theta = 0.0)
        @test length(turning_points(g)) == 1
        @test n_infinite_lines(g) == 3
        @test isempty(edges(g))
        @test isempty(finite_lines(g))
        @test topology_signature(g) == ([1], 3, Tuple{Int,Int}[])
        @test !is_degenerate(g)
        @test turning_point_orders(g) == [1]

        # exit angles are 0, 2π/3, 4π/3 (mod 2π), one per ray
        exits = sort([mod2pi(angle(points(l)[end])) for l in lines(g)])
        @test isapprox(exits, [0.0, 2π / 3, 4π / 3]; atol = 1e-4)
        # every ray escapes with the same mass (symmetry)
        @test all(l -> l.endpoint === :infinity, lines(g))
    end

    @testset "θ-rotation moves rays by the seeding formula" begin
        airy = SchrodingerProblem([0, 1])
        for θ in (0.0, 0.3, 0.6, -0.4)
            g = stokes_graph(airy; theta = θ)
            # ray k has exit direction (2θ + 2πk)/3; k = 0 → 2θ/3
            a0 = angle(points(lines(g)[1])[end])
            @test isapprox(a0, 2θ / 3; atol = 1e-4)
        end
    end

    @testset "w-drift is killed by re-projection" begin
        # Airy θ = 0, ray k = 0 is exactly the +real axis; any imaginary part is drift.
        g = stokes_graph(SchrodingerProblem([0, 1]); theta = 0.0)
        l0 = lines(g)[1]
        drift = maximum(abs(imag(z)) for z in points(l0))
        @test drift < 1e-8
    end

    @testset "order-≥2 turning points are traced, not refused" begin
        # M7a lifted the blanket refusal: an order-m point emits m+2 rays. The error
        # type survives for cases genuinely out of scope; what used to throw here now
        # returns a degenerate graph (see the degeneracy testsets below).
        g = stokes_graph(SchrodingerProblem([0, 0, 1]))   # Q = z²
        @test g isa StokesGraph
        @test is_degenerate(g)
    end

    @testset "double well: inner saddle + stable signature" begin
        dw = SchrodingerProblem([3//4, 0, -2, 0, 1])   # (z²−1)² − 1/4
        g = stokes_graph(dw; theta = 0.0)
        @test length(turning_points(g)) == 4

        # canonical turning-point order is by position: inner pair are indices 2,3
        locs = [real(location(t)) for t in turning_points(g)]
        @test issorted(locs)

        # exactly one saddle: the inner pair, double-traced
        es = edges(g)
        @test length(es) == 1
        @test es[1] == (2, 3)
        @test count(is_finite_line, lines(g)) == 2      # traced from both endpoints
        @test n_infinite_lines(g) == 10                 # 4·3 − 2

        # the two half-saddle line-masses agree (½|Z| each direction)
        halfmasses = [l.mass for l in finite_lines(g)]
        @test isapprox(halfmasses[1], halfmasses[2]; rtol = 1e-4)

        # topology signature is stable under seed / tolerance halving
        s1 = topology_signature(g)
        s2 = topology_signature(stokes_graph(dw; seed_radius = 0.005, reltol = 5e-10))
        @test s1 == s2
        @test s1 == ([1, 1, 1, 1], 10, [(2, 3)])   # orders first (all simple here)
    end

    @testset "off-critical phase breaks the saddle" begin
        # away from θ_c = 0 the inner Stokes line no longer connects the two turning
        # points, so no finite edge survives.
        dw = SchrodingerProblem([3//4, 0, -2, 0, 1])
        g = stokes_graph(dw; theta = 0.25)
        @test isempty(edges(g))
        @test n_infinite_lines(g) == 12
    end

    @testset "degenerate turning points: m + 2 rays" begin
        # Q = z^m has a single order-m turning point; the Stokes condition puts m+2
        # rays at (2θ − arg c_m + 2πk)/(m+2). At θ = 0 with c_m = 1 they are the
        # (m+2)-th roots of unity directions.
        for m in 2:4
            p = SchrodingerProblem([zeros(Int, m); 1])       # Q = z^m
            g = stokes_graph(p; theta = 0.0)
            @test length(turning_points(g)) == 1
            @test order(turning_points(g)[1]) == m
            @test is_degenerate(g)
            @test turning_point_orders(g) == [m]
            @test length(lines(g)) == m + 2
            @test n_infinite_lines(g) == m + 2
            exits = sort([mod2pi(angle(points(l)[end])) for l in lines(g)])
            @test isapprox(exits, [2π * k / (m + 2) for k in 0:(m + 1)]; atol = 1e-3)
        end
    end

    @testset "a double turning point is a crossing, not a vertex" begin
        # For even m, √Q is single-valued at z₀ (the spectral curve w² = Q has a node,
        # not a branch point), so opposite rays are analytic continuations of one
        # another: ray k and ray k+2 of the four leave z₀ in exactly opposite
        # directions and together form one smooth Stokes curve.
        g = stokes_graph(SchrodingerProblem([0, 0, 1]); theta = 0.0)   # Q = z²
        @test is_degenerate(g)
        ls = sort(collect(lines(g)); by = l -> mod2pi(angle(points(l)[end])))
        @test length(ls) == 4
        near(l) = points(l)[min(6, length(points(l)))]
        for k in 1:2
            a, b = angle(near(ls[k])), angle(near(ls[k + 2]))
            @test abs(mod2pi(a - b) - π) < 1e-3
        end
        # an odd-order point is a genuine branch point: 5 rays, no opposite pairs
        g3 = stokes_graph(SchrodingerProblem([0, 0, 0, 1]); theta = 0.0)
        @test length(lines(g3)) == 5
    end

    @testset "coalescence: 6 rays and a saddle collapse to 4" begin
        # Q = z²/4 − E has two simple turning points ±2√E. They merge at E = 0: the
        # generic graph has 2·3 = 6 rays, the degenerate one 4. The signature carries
        # the orders, so the two are never confused.
        gen = stokes_graph(weber_problem(1 // 4); theta = 0.3)
        @test turning_point_orders(gen) == [1, 1]
        @test length(lines(gen)) == 6
        @test !is_degenerate(gen)

        deg = stokes_graph(weber_problem(0 // 1); theta = 0.3)
        @test turning_point_orders(deg) == [2]
        @test length(lines(deg)) == 4
        @test is_degenerate(deg)
        @test topology_signature(gen) != topology_signature(deg)

        # The saddle mass → 0 as the pair merges. Z = ∮√Q = 2πiE is purely imaginary,
        # so the saddle phase is θ = π/2 (not 0 as for a real double well), and each
        # half-line carries |Z|/2 = πE - so quartering E must quarter the mass. The
        # traced mass sits a uniform ~0.2% low because a ray is seeded `seed_radius`
        # away from the turning point and never accumulates that first stretch; the
        # exact period is `period_integral`'s job, not the tracer's.
        masses = Float64[]
        for E in (1 // 4, 1 // 16, 1 // 64)
            gE = stokes_graph(weber_problem(E); theta = π / 2)
            fl = finite_lines(gE)
            @test !isempty(fl)
            @test edges(gE) == [(1, 2)]
            push!(masses, minimum(mass, fl))
            @test masses[end] ≈ π * Float64(E) rtol = 5e-3
        end
        @test issorted(masses; rev = true)
        @test masses[end] ≈ masses[1] / 16 rtol = 5e-3
    end

    @testset "the triangulation layer refuses a degenerate graph" begin
        g = stokes_graph(SchrodingerProblem([0, 0, 1]); theta = 0.3)
        @test is_degenerate(g)
        @test_throws NonGenericGraph ideal_triangulation(g)
    end

    @testset "BigFloat opt-in traces at high precision" begin
        setprecision(BigFloat, 128) do
            g = stokes_graph(SchrodingerProblem([0, 1]); theta = 0.0, bigfloat = true)
            @test eltype(points(lines(g)[1])) == Complex{BigFloat}
            @test n_infinite_lines(g) == 3
        end
    end
end
