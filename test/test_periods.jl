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

    @testset "a pole must not inflate the turning-point tolerance" begin
        # Regression. `_period_atol` scaled `sqrt(eps)` by the MAXIMUM |Q| on the path.
        # That is not a typical value once Q has poles: a contour running near one sent
        # the tolerance far above honest |Q| values, and a perfectly good vertex was
        # rejected as a turning point - which broke `charge_basis(prob, g)` on two of
        # three chambers of this very problem.
        p = RationalProblem([-0.3, -1.0, 0.0, 1.0], [0.0], [4])
        g = stokes_graph(p; theta = 0.9)
        c = charge_contour(p, diagonal_core_paths(g)[1])
        # the contour really does dip next to the order-4 pole, so max|Q| really is huge
        @test minimum(abs, c) < 0.05
        @test maximum(abs(p(z)) for z in c) > 1e5
        # the invariant that was violated: the tolerance must sit below every |Q| the
        # contour actually takes, or it flags vertices that are nowhere near a zero
        @test ExactWKB._period_atol(p, c, Float64) < minimum(abs(p(z)) for z in c)
        @test isfinite(period_integral(p, c))
    end

    @testset "parameter derivatives (frozen contour)" begin
        @testset "harmonic ∂/∂E ∮√(z²−E) = −iπ" begin
            for E in (1.0, 2.5, 4.0)
                harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = E)
                tps = turning_points(harm)
                c = encircling_contour(tps[1], tps[2])
                @test period_derivative(harm, c) ≈ -im * π rtol = 1e-8
                w = wkb_expansion(harm; order = 3)
                @test wkb_period(wkb_derivative(w), c, -1) ≈ -im * π rtol = 1e-8
            end
        end

        @testset "the two routes agree (Riccati ∂S₋₁ vs ∂Q/2√Q)" begin
            # `period_derivative` integrates ∂Q/(2√Q) directly; `wkb_period` of a
            # WKBDerivative integrates the term algebra's ∂S₋₁. Same object, two paths.
            prob = SchrodingerProblem([0 // 1, 0 // 1, -4 // 1, 0 // 1, 1 // 1];
                                      energy = -2 // 1)
            c = spectral_cycles(prob, -2 // 1)[1].contour
            w = wkb_expansion(prob; order = 3)
            for (wrt, idx) in ((:energy, 0), (:coefficient, 0), (:coefficient, 1),
                               (:coefficient, 3), (:coefficient, 4))
                a = period_derivative(prob, c; wrt, index = idx)
                b = wkb_period(wkb_derivative(w; wrt, index = idx), c, -1)
                @test a ≈ b rtol = 1e-20
            end
        end

        @testset "quantum orders against a BigFloat finite difference" begin
            # The FD is taken at the SAME frozen contour, so this isolates the term
            # algebra from the contour question that the ledger settles.
            setprecision(BigFloat, 128) do
                v = [0 // 1, 0 // 1, -4 // 1, 0 // 1, 1 // 1]
                E0 = -2 // 1
                prob = SchrodingerProblem(v; energy = E0)
                c = spectral_cycles(prob, E0)[1].contour
                w = wkb_expansion(prob; order = 5)
                dw = wkb_derivative(w)
                h = big(1) / big(10)^10
                wp = wkb_expansion(with_energy(prob, E0 + h); order = 5)
                wm = wkb_expansion(with_energy(prob, E0 - h); order = 5)
                for m in (-1, 1, 3, 5)
                    ex = wkb_period(dw, c, m)
                    fd = (wkb_period(wp, c, m) - wkb_period(wm, c, m)) / (2h)
                    @test abs(ex - fd) < 1e-18 * (1 + abs(ex))
                end
            end
        end

        @testset "∂/∂v_k against a coefficient finite difference" begin
            setprecision(BigFloat, 128) do
                v = [0 // 1, 0 // 1, -4 // 1, 0 // 1, 1 // 1]
                E0 = -2 // 1
                prob = SchrodingerProblem(v; energy = E0)
                c = spectral_cycles(prob, E0)[1].contour
                h = big(1) / big(10)^10
                for k in 0:4
                    ex = period_derivative(prob, c; wrt = :coefficient, index = k)
                    vp = collect(v); vp[k + 1] += h
                    vm = collect(v); vm[k + 1] -= h
                    fd = (period_integral(SchrodingerProblem(vp; energy = E0), c) -
                          period_integral(SchrodingerProblem(vm; energy = E0), c)) / (2h)
                    @test abs(ex - fd) < 1e-16 * (1 + abs(ex))
                end
            end
        end

        @testset "ledger item 2: refused on an open contour" begin
            harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = 1.0)
            tps = turning_points(harm)
            c = encircling_contour(tps[1], tps[2])
            w = wkb_expansion(harm; order = 2)
            @test_throws ContourError period_derivative(harm, c; closed = false)
            @test_throws ContourError wkb_period(wkb_derivative(w), c, -1; closed = false)
            # the primal is unaffected: an open contour is still a legal period
            @test isfinite(period_integral(harm, c; closed = false))
            @test isfinite(wkb_period(w, c, -1; closed = false))
        end

        @testset "argument checks" begin
            harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = 1.0)
            tps = turning_points(harm)
            c = encircling_contour(tps[1], tps[2])
            @test_throws Resurgence.InvalidArgument period_derivative(harm, c;
                                                                      wrt = :nonsense)
            @test_throws Resurgence.InvalidArgument period_derivative(
                harm, c; wrt = :coefficient, index = 3)   # V has degree 2
        end
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
