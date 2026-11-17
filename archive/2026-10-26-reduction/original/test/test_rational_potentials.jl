# Rational potentials: Q with poles. The oracles are (a) the algebra of Q = N/D
# against a hand-differentiated closed form, (b) the local structure at a pole - the
# asymptotic directions and the Stokes rays that run into it - and (c) the payoff, the
# Mathieu / pure-SU(2) problem, whose turning points, surface and periods are all
# known independently.
#
# The last testset is the four-package loop: a Mathieu potential's Stokes graph is
# traced on the annulus, dualised to the Kronecker quiver, and its two classical
# periods are compared with `sw_curve.jl`'s `a` and `a_D` - the SAME theory computed
# by two disjoint routes (Schrödinger tracer + branch-tracked quadrature vs
# Seiberg-Witten curve + Carlson elliptic integrals), which had never been compared.

import ClusterAlgebras

@testset "rational potentials" begin
    @testset "construction and validation" begin
        # Q = (1 + z²)/z³
        p = RationalProblem([1.0, 0.0, 1.0], [0.0], [3])
        @test q_numerator(p) == [1.0, 0.0, 1.0]
        @test poles(p) == [0.0]
        @test pole_orders(p) == [3]
        @test n_finite_poles(p) == 1
        @test variable(p) == :z
        @test degree(p) == 2 - 3            # the exponent of Q at infinity
        @test p == RationalProblem([1.0, 0.0, 1.0], [0.0], [3])
        @test p != RationalProblem([1.0, 0.0, 2.0], [0.0], [3])

        # a polynomial problem has no poles, and `degree` means the same thing
        sp = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
        @test isempty(poles(sp))
        @test isempty(pole_orders(sp))
        @test n_finite_poles(sp) == 0
        @test sp isa AbstractSchrodingerProblem
        @test p isa AbstractSchrodingerProblem

        # a double pole is a puncture and IS supported; a simple pole is an orbifold
        # point and is not, and infinity must stay irregular
        @test RationalProblem([1.0, 1.0], [0.0], [2]) isa RationalProblem
        @test_throws InvalidPotential RationalProblem([1.0, 1.0], [0.0], [1])
        # Q ~ z^{-2} at ∞ ⇒ pole order 2 there
        @test_throws InvalidPotential RationalProblem([1.0], [0.0], [3])
        @test_throws InvalidPotential RationalProblem(Float64[], [0.0], [3])
        @test_throws InvalidPotential RationalProblem([0.0, 0.0], [0.0], [3])
        @test_throws InvalidPotential RationalProblem([1.0, 1.0, 1.0], [0.0, 0.0], [3, 3])
        @test_throws InvalidPotential RationalProblem([1.0, 1.0, 1.0], [0.0], [3, 3])
        # the numerator vanishing at a pole means the declared order is wrong
        @test_throws InvalidPotential RationalProblem([0.0, 1.0, 1.0, 1.0], [0.0], [3])
    end

    @testset "evaluation, derivative, Taylor shift" begin
        # Q = (1 + z²)/z³ = z^{-3} + z^{-1};  Q′ = −3z^{-4} − z^{-2}
        p = RationalProblem([1.0, 0.0, 1.0], [0.0], [3])
        for z in (0.7, -1.3, 2.0 + 0.5im)
            @test p(z) ≈ z^-3 + z^-1
            @test q_derivative_at(p, z) ≈ -3z^-4 - z^-2
        end

        # two poles, so the denominator is a genuine product (numerator degree 6 keeps
        # infinity irregular: e = 6 − 7 = −1)
        p2 = RationalProblem([1.0, 2.0, 3.0, 4.0, 1.0, 0.0, 2.0], [1.0, -2.0], [3, 4])
        N(z) = 1 + 2z + 3z^2 + 4z^3 + z^4 + 2z^6
        D(z) = (z - 1)^3 * (z + 2)^4
        for z in (0.3, 2.5, -0.4 + 1.1im)
            @test p2(z) ≈ N(z) / D(z)
            # central difference against the analytic quotient rule
            h = 1e-6
            @test q_derivative_at(p2, z) ≈ (p2(z + h) - p2(z - h)) / 2h rtol = 1e-6
        end

        # the Taylor shift reproduces Q on a small disk about z0
        z0 = 0.4
        tay = q_taylor_at(p2, z0)
        for δ in (0.01, -0.02, 0.015im)
            @test sum(tay[k] * δ^(k - 1) for k in eachindex(tay)) ≈ p2(z0 + δ) rtol = 1e-6
        end
        @test tay[1] ≈ p2(z0)
        @test tay[2] ≈ q_derivative_at(p2, z0)
        @test_throws InvalidPotential q_taylor_at(p2, 1.0)     # z0 is a pole
    end

    @testset "turning points are the numerator's zeros" begin
        # Q = (z² − 1)/z³: turning points ±1, and the poles are NOT turning points
        p = RationalProblem([-1.0, 0.0, 1.0], [0.0], [3])
        tps = turning_points(p)
        @test length(tps) == 2
        @test sort(real.(location.(tps))) ≈ [-1.0, 1.0]
        @test all(is_simple, tps)
        @test all(t -> abs(p(location(t))) < 1e-10, tps)
    end

    @testset "asymptotic directions" begin
        # `Q ~ c ζ^e` puts |e + 2| directions at (2θ − arg c + 2πk)/(e + 2)
        p = RationalProblem([1.0, 0.0, 1.0], [0.0], [3])       # e = −1 at ∞, m = 3
        for θ in (0.0, 0.3, 1.1)
            d∞ = asymptotic_directions(p, θ)
            dp = asymptotic_directions(p, θ; pole = 1)
            @test length(d∞) == 1                              # e + 2 = 1
            @test length(dp) == 1                              # 2 − m = −1
            @test d∞[1] ≈ mod(2θ, 2π)                          # arg c = 0
            @test dp[1] ≈ mod(-2θ, 2π)
        end
        # a polynomial problem: the familiar d + 2 asymptotic directions
        cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
        @test length(asymptotic_directions(cubic, 0.3)) == 5
        # a higher-order pole carries m − 2 marked points
        p5 = RationalProblem([1.0, 0.0, 0.0, 0.0, 1.0], [0.0], [5])   # e = 4 − 5 = −1
        @test length(asymptotic_directions(p5, 0.2; pole = 1)) == 3
        @test length(asymptotic_directions(p5, 0.2)) == 1
    end

    @testset "the tracer runs rays into poles" begin
        p = mathieu_problem(1.0, 3.0)
        g = stokes_graph(p; theta = 0.3)
        @test length(ExactWKB.lines(g)) == 6                   # 3 rays × 2 turning points
        @test isempty(ExactWKB.finite_lines(g))                # θ = 0.3 is off the wall
        @test length(ExactWKB.boundary_lines(g)) == 6
        @test n_infinite_lines(g) == 3
        @test count(l -> ExactWKB.endpoint(l) === :pole, ExactWKB.lines(g)) == 3
        @test all(l -> ExactWKB.target(l) == 1,
                  filter(l -> ExactWKB.endpoint(l) === :pole, ExactWKB.lines(g)))
        @test ExactWKB.poles(g) == ComplexF64[0.0]
        @test ExactWKB.pole_orders(g) == [3]

        # every ray lands near its circle's single asymptotic direction
        φ∞ = asymptotic_directions(p, 0.3)[1]
        φp = asymptotic_directions(p, 0.3; pole = 1)[1]
        for l in ExactWKB.lines(g)
            z = ExactWKB.points(l)[end]
            if ExactWKB.endpoint(l) === :infinity
                @test abs(ExactWKB._wrap(angle(z) - φ∞)) < π
            else
                @test abs(ExactWKB._wrap(angle(z) - φp)) < π
            end
        end
    end

    @testset "mathieu_problem is the SU(2) curve on the w-plane" begin
        for (Λ, u) in ((1.0, 3.0), (1.2, 5.0), (0.7, 1.5))
            p = mathieu_problem(Λ, u)
            @test n_finite_poles(p) == 1
            @test pole_orders(p) == [3]
            @test degree(p) == -1                              # pole order 3 at ∞ too
            # the two turning points in closed form; they collide at u = ±2Λ²
            r = sqrt(u^2 - 4Λ^4)
            w = sort([(u - r) / (2Λ^2), (u + r) / (2Λ^2)])
            @test sort(real.(location.(turning_points(p)))) ≈ w
            @test prod(w) ≈ 1                                  # w₁w₂ = 1
            # Q̃(w) = −(Λ²w² − u w + Λ²)/(4w³)
            for z in (0.6, 1.7, 0.3 + 0.2im)
                @test p(z) ≈ -(Λ^2 * z^2 - u * z + Λ^2) / (4z^3)
            end
        end
        # at the Seiberg-Witten singularity the turning points coalesce
        pc = mathieu_problem(1.0, 2.0)
        @test length(turning_points(pc)) == 1
        @test order(turning_points(pc)[1]) == 2
    end

    @testset "the four-package loop: WKB side vs curve side of SU(2)" begin
        # (1) the Stokes graph of the Mathieu potential triangulates the ANNULUS
        p = mathieu_problem(1.0, 3.0)
        t = ideal_triangulation(stokes_graph(p; theta = 0.3))
        @test ExactWKB.n_boundaries(t) == 2
        @test n_marked_points(t) == 2
        @test t.marked_boundary == [1, 2]                   # one point per circle
        @test n_diagonals(t) == 2
        # both arcs join the two marked points; the two boundary segments are loops
        @test t.edge_endpoints == [(1, 2), (1, 2), (1, 1), (2, 2)]
        @test t.is_diagonal == [true, true, false, false]

        # (2) the quiver is Ã(1,1) = Kronecker = the SU(2) BPS quiver
        B = triangulation_quiver(t).B
        @test abs.(B) == [0 2; 2 0]
        @test ClusterAlgebras.canonical_form(ClusterAlgebras.Quiver(B)).B ==
              ClusterAlgebras.canonical_form(su2_bps_quiver()).B
        @test string(ClusterAlgebras.mutation_type(ClusterAlgebras.Quiver(B))) ==
              "A(1,1)^(1)"

        # (3) flip = mutation, on the annulus as on the disk
        for k in 1:2
            f = flip(t, k)
            @test triangulation_quiver(f).B ==
                  ClusterAlgebras.mutate(ClusterAlgebras.Quiver(B), k).B
            back = flip(f, k; direction = -1)
            @test back.edge_endpoints == t.edge_endpoints
            @test back.triangles == t.triangles
            @test back.diagonal_tp_pair == t.diagonal_tp_pair
        end

        # (4) the central charges. Under w = e^{2ix} the electric cycle is |w| = 1 and
        # the magnetic cycle encircles the two turning points; both periods are
        # ∓(1/iπ)∮√Q̃ dw (see the ledger in src/rational_potentials.jl). Nothing here
        # touches sw_curve.jl's machinery - it is elliptic integrals on the SW curve,
        # this is a branch-tracked quadrature on the Schrödinger problem.
        circle = [cis(2π * j / 512) for j in 0:511]
        for (Λ, u) in ((1.0, 3.0), (1.2, 5.0), (1.0, 10.0), (0.7, 1.5))
            prob = mathieu_problem(Λ, u)
            zs = location.(turning_points(prob))
            a = period_integral(prob, circle) / (im * π)
            m = 0.4 * min(abs(zs[1]), abs(zs[2] - zs[1]) / 2)
            aD = -period_integral(prob, encircling_contour(zs[1], zs[2];
                                                           margin = m, n = 256)) / (im * π)
            per = sw_periods(SeibergWittenSU2(Λ = Λ), u)
            @test a ≈ per.a atol = 1e-12
            @test aD ≈ per.a_D atol = 1e-12
        end

        # (5) the saddle the tracer finds is the magnetic (monopole) state: at u = 3,
        # Λ = 1 the wall sits at θ = 0 and |Z| is the magnetic period π·|a_D|. The
        # tolerance is the straight-line saddle path of `saddle_candidates`, not the
        # period layer - the contour integral in (4) agrees to 1e-12.
        sads = saddle_candidates(p)
        @test length(sads) == 1
        @test ExactWKB.mass(sads[1]) ≈
              π * abs(sw_periods(SeibergWittenSU2(Λ = 1.0), 3.0).a_D) rtol = 1e-7
        @test abs(ExactWKB.theta(sads[1])) < 1e-9

        # (6) the u-DERIVATIVES of the same two periods, differentiated under the
        # integral sign at frozen contour. ∂_u Q̃ = 1/(4w²), and a closed cycle's period
        # does not depend on its contour, so nothing has to track the moving turning
        # points. Three independent routes have to land on one number: the closed-form
        # elliptic K (`sw_period_derivatives`), the Picard-Fuchs transport
        # (`continue_periods`), and this quadrature.
        dQdu(w) = 1 / (4 * w^2)
        for (Λ, u) in ((1.0, 3.0), (1.2, 5.0), (1.0, 10.0), (0.7, 1.5))
            prob = mathieu_problem(Λ, u)
            zs = location.(turning_points(prob))
            m = 0.4 * min(abs(zs[1]), abs(zs[2] - zs[1]) / 2)
            ct = encircling_contour(zs[1], zs[2]; margin = m, n = 256)
            da = period_derivative(prob, circle, dQdu) / (im * π)
            daD = -period_derivative(prob, ct, dQdu) / (im * π)
            ref = sw_period_derivatives(SeibergWittenSU2(Λ = Λ), u)
            @test da ≈ ref.da atol = 1e-12
            @test daD ≈ ref.da_D atol = 1e-12
        end

        # the Picard-Fuchs route, transported along a path that leaves the real axis
        cp = continue_periods(SeibergWittenSU2(Λ = 1.0),
                              [3.0 + 0.0im, 4.0 + 0.5im, 5.0 + 0.0im])
        p5 = mathieu_problem(1.0, 5.0)
        z5 = location.(turning_points(p5))
        ct5 = encircling_contour(z5[1], z5[2];
                                 margin = 0.4 * min(abs(z5[1]), abs(z5[2] - z5[1]) / 2),
                                 n = 256)
        @test cp.da ≈ period_derivative(p5, circle, dQdu) / (im * π) atol = 1e-10
        @test cp.da_D ≈ -period_derivative(p5, ct5, dQdu) / (im * π) atol = 1e-10

        # Cauchy-Riemann: Z is holomorphic in u, so the real and imaginary difference
        # quotients must give the same derivative. Taken on the electric cycle, whose
        # first vertex w = 1 has Q̃ = (u − 2Λ²)/4 > 0 - the magnetic contour starts where
        # arg Q̃ ≈ π, and a finite difference there flips the √Q branch. That fragility
        # is the finite difference's, not the derivative's.
        let Λ = 1.0, u = 3.0, h = 1e-5
            ex = period_derivative(mathieu_problem(Λ, u), circle, dQdu)
            fre = (period_integral(mathieu_problem(Λ, u + h), circle) -
                   period_integral(mathieu_problem(Λ, u - h), circle)) / (2h)
            fim = (period_integral(mathieu_problem(Λ, u + im * h), circle) -
                   period_integral(mathieu_problem(Λ, u - im * h), circle)) / (2im * h)
            @test ex ≈ fre atol = 1e-9
            @test ex ≈ fim atol = 1e-9
        end
    end

    @testset "the charge lattice on the annulus" begin
        # A turning-point pair is not a name for a cycle here: BOTH diagonals join
        # turning points 1 and 2, and they differ by the loop around the hole. So
        # `charge_basis(prob, t)` cannot express them and `charge_basis(prob, g)` -
        # which builds each cycle from its own strip region - must.
        Λ, u = 1.0, 3.0
        p = mathieu_problem(Λ, u)
        g = stokes_graph(p; theta = 0.3)
        t = ideal_triangulation(g)
        @test t.diagonal_tp_pair == [(1, 2), (1, 2)]        # the ambiguity, explicitly

        paths = diagonal_core_paths(g)
        @test length(paths) == 2
        zs = location.(turning_points(p))
        for path in paths
            @test path[1] ≈ zs[1]
            @test path[end] ≈ zs[2]
        end
        # the two strips are genuinely different classes: one wraps, one does not
        @test maximum(abs, paths[1]) > 1 && maximum(abs, paths[2]) > 1
        @test minimum(real, paths[1]) != minimum(real, paths[2])

        cb = charge_basis(p, g)
        @test n_charges(cb) == 2
        # the keystone, which is what actually pins the homotopy classes - winding
        # numbers are homology downstairs and cannot see the wrap
        @test signed_pairing(cb) == -triangulation_quiver(t).B
        @test abs.(cb.pairing) == [0 2; 2 0]

        # the physical content: the two basis cycles are ±monopole and ±dyon of the
        # SU(2) lattice, with Z_contour = iπ·(n_m a_D + n_e a) from sw_curve.jl
        per = sw_periods(SeibergWittenSU2(Λ = Λ), u)
        Zmon = im * π * per.a_D                             # (1, 0)
        Zdyon = im * π * (-per.a_D + 2 * per.a)             # (−1, 2)
        Z = physical_charges(cb)
        @test minimum(abs, [abs(Z[1] - Zmon), abs(Z[1] + Zmon),
                            abs(Z[1] - Zdyon), abs(Z[1] + Zdyon)]) < 1e-11
        @test minimum(abs, [abs(Z[2] - Zmon), abs(Z[2] + Zmon),
                            abs(Z[2] - Zdyon), abs(Z[2] + Zdyon)]) < 1e-11
        # one of each, not the same state twice
        @test (abs(abs(Z[1]) - abs(Zmon)) < 1e-11) != (abs(abs(Z[2]) - abs(Zmon)) < 1e-11)

        # the bridge seed is the Kronecker exchange matrix
        @test bridge_seed(cb).quiver.B == triangulation_quiver(t).B
    end

    @testset "a second punctured surface: Ã(2,1)" begin
        # One order-4 pole plus infinity ⇒ 2 marked points on one boundary circle and 1
        # on the other: the annulus Ã(2,1), rank 3. Everything above only ever exercised
        # Ã(1,1), where both circles carry a single marked point and every strip has a
        # distinct side - too special. This fixture is what found the `_period_atol`
        # defect (see test_periods.jl), so it is here to stay.
        p = RationalProblem([-0.3, -1.0, 0.0, 1.0], [0.0], [4])
        @test pole_orders(p) == [4]
        @test degree(p) == -1
        @test length(asymptotic_directions(p, 0.4)) == 1            # e + 2 at infinity
        @test length(asymptotic_directions(p, 0.4; pole = 1)) == 2  # m − 2 at the pole
        @test length(turning_points(p)) == 3
        @test all(is_simple, turning_points(p))

        for θ in (0.4, 0.9, 1.7)
            g = stokes_graph(p; theta = θ)
            t = ideal_triangulation(g)
            @test n_boundaries(t) == 2
            @test n_marked_points(t) == 3
            @test sort(t.marked_boundary) == [1, 2, 2]
            @test n_diagonals(t) == 3                    # (3n − M)/2 = (9 − 3)/2
            @test length(t.triangles) == 3
            B = triangulation_quiver(t).B
            @test string(ClusterAlgebras.mutation_type(ClusterAlgebras.Quiver(B))) ==
                  "A(1,2)^(1)"
            # the charge lattice closes in every chamber, keystone and all
            cb = charge_basis(p, g)
            @test n_charges(cb) == 3
            @test signed_pairing(cb) == -B
            @test all(!iszero, physical_charges(cb))
        end
    end
end
