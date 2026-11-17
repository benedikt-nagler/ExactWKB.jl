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
        p = RationalProblem([-0.25, 0.75, -0.25], [0.0], [3])
        g = stokes_graph(p; theta = 0.3)
        @test length(ExactWKB.lines(g)) == 6                   # 3 rays × 2 turning points
        @test isempty(ExactWKB.finite_lines(g))                # θ = 0.3 is off the wall
        @test count(l -> endpoint(l) in (:pole, :infinity), ExactWKB.lines(g)) == 6
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


end
