using QuadGK: quadgk

@testset "periods" begin
    @testset "harmonic ∮√(z²−E) = −iπE (orientation/branch)" begin
        for E in (1.0, 2.5, 4.0)
            harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = E)
            tps = turning_points(harm)
            P = period_integral(harm, encircling_contour(tps[1], tps[2]))
            @test P ≈ -im * π * E rtol = 1e-8
        end
    end

    @testset "elliptic A-cycle vs independent real integral" begin
        # Q = (z²−1)(z²−k²); cycle around the pair (k, 1). On (k,1) Q<0, so
        # |∮√Q| = 2∫_k^1 √|Q| dx - an independent real integral.
        k = 0.6
        # Q = z⁴ − (1+k²)z² + k²
        prob = SchrodingerProblem([k^2, 0.0, -(1 + k^2), 0.0, 1.0])
        tps = sort(turning_points(prob); by = t -> real(location(t)))
        # turning points at −1, −k, k, 1 → the (k,1) pair are tps[3], tps[4]
        P = period_integral(prob, encircling_contour(tps[3], tps[4]))
        ref, _ = quadgk(x -> sqrt(abs((x^2 - 1) * (x^2 - k^2))), k, 1.0)
        @test abs(P) ≈ 2 * ref rtol = 1e-6
        @test abs(real(P)) < 1e-6 * abs(P)          # purely imaginary
    end

    @testset "wkb_period reproduces the classical period at m = -1" begin
        harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = 1.0)
        w = wkb_expansion(harm; order = 4)
        tps = turning_points(harm)
        c = encircling_contour(tps[1], tps[2])
        @test wkb_period(w, c, -1) ≈ period_integral(harm, c) rtol = 1e-8
        @test_throws Resurgence.InvalidArgument wkb_period(w, c, 5)
    end

    @testset "ContourError cases" begin
        harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = 1.0)
        tps = turning_points(harm)
        # a vertex sitting on a turning point
        bad = [tps[1].z, 1.0 + 1.0im, -1.0 + 1.0im]
        @test_throws ContourError period_integral(harm, bad; closed = true)
        # a loop around a single turning point → √Q not single-valued
        loop = [tps[2].z + 0.3 * cis(2π * j / 24) for j in 0:23]
        @test_throws ContourError period_integral(harm, loop; closed = true)
    end

    @testset "BigFloat contour → BigFloat result" begin
        setprecision(256) do
            E = big"2.0"
            harm = SchrodingerProblem([big(0.0), big(0.0), big(1.0)]; energy = E)
            tps = turning_points(harm)
            c = encircling_contour(tps[1], tps[2])
            @test eltype(c) == Complex{BigFloat}
            P = period_integral(harm, c)
            @test P isa Complex{BigFloat}
            @test abs(P - (-im * BigFloat(π) * E)) < 1e-25
        end
    end
end
