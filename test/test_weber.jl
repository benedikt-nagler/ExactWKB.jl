using Resurgence: FormalSeries, coefficients, power_offset, borel, pade, poles

@testset "weber" begin
    @testset "Bernoulli numbers" begin
        B = ExactWKB._bernoulli(12)
        @test B[1] == 1
        @test B[2] == -1 // 2
        @test B[3] == 1 // 6            # B₂
        @test B[5] == -1 // 30          # B₄
        @test B[7] == 1 // 42           # B₆
        @test B[9] == -1 // 30          # B₈
        @test B[11] == 5 // 66          # B₁₀
        @test B[13] == -691 // 2730     # B₁₂
        @test all(iszero, B[4:2:13])    # odd-index B vanish beyond B₁
    end

    @testset "exact Voros series = the log Γ(½+ν) Stirling correction" begin
        S = weber_voros_series(; order = 9)
        c = coefficients(S)
        @test power_offset(S) == 1
        @test eltype(c) <: Rational
        # c_n = (2^{1−2n} − 1) B_{2n} / (2n(2n−1)); DLMF §5.11
        @test c[1] == -1 // 24
        @test c[3] == 7 // 2880
        @test c[5] == -31 // 40320
        @test c[7] == 127 // 215040
        @test all(iszero, c[2:2:end])   # odd-only, like every Voros series here
        @test_throws Resurgence.InvalidArgument weber_voros_series(; order = 0)
    end

    @testset "weber_model: index, orientation, refusals" begin
        # Normal form Q = z²/4 − E: turning points ±2√E, and ν = E/ħ exactly.
        for E in (1 // 1, 4 // 1, 1 // 4)
            p = weber_problem(E)
            m = weber_model(p)
            a = 2 * sqrt(Float64(E))
            @test sort(Float64.(real.(collect(m.turning_points)))) ≈ [-a, a] rtol = 1e-20
            # ledger item 1: Z = 2πiE and ν = E/ħ (positive on the real-E branch)
            @test Complex{Float64}(merging_period(m)) ≈ 2im * π * Float64(E) rtol = 1e-12
            for ħ in (1.0, 0.25)
                @test Complex{Float64}(weber_index(m, ħ)) ≈ Float64(E) / ħ rtol = 1e-12
            end
        end
        # a general merging pair (the quartic double well) is accepted
        dw = SchrodingerProblem([3 // 4, 0, -2, 0, 1])
        @test weber_model(dw) isa WeberModel
        # already merged: ν is not defined there
        @test_throws CoalescentTurningPoints weber_model(SchrodingerProblem([0, 0, 1]))
    end

    @testset "computed Voros coefficient reproduces the closed form" begin
        # THE headline oracle of this layer: the package's own WKB recursion,
        # integrated numerically out to infinity, is the Stirling series for
        # log Γ(½ + ν) (Takei 2008, Koike–Takei 2011).
        p = weber_problem(1 // 1)
        m = weber_model(p)
        pred = coefficients(weber_voros_series(m; order = 9))
        got = coefficients(weber_voros_coefficient(p, m; order = 9))
        @test length(got) == length(pred)
        @test all(iszero, got[2:2:end])
        for k in 1:2:9
            @test Complex{Float64}(got[k]) ≈ Complex{Float64}(pred[k]) rtol = 1e-6
        end
        # sharper on the tail orders, where the R-extrapolation is not the floor
        for k in 3:2:9
            @test Complex{Float64}(got[k]) ≈ Complex{Float64}(pred[k]) rtol = 1e-9
        end

        # Ledger item 2 (the factor 2) is a normalization, not an error: dropping it
        # would miss uniformly at EVERY order, so halving must break all of them.
        @test all(k -> !isapprox(Complex{Float64}(got[k]) / 2,
                                 Complex{Float64}(pred[k]); rtol = 1e-3), 1:2:9)

        # Ledger item 1 (the sign of ν): W has odd powers of 1/ν, so ν → −ν flips
        # every coefficient - the comparison fails on the whole vector.
        flipped = coefficients(weber_voros_series(
            WeberModel((m.turning_points), -merging_period(m), -m.index_scale); order = 9))
        @test all(k -> !isapprox(Complex{Float64}(got[k]), Complex{Float64}(flipped[k]);
                                 rtol = 1e-3), 1:2:9)

        # Ledger item 5: the other turning point of the pair is a genuinely different
        # ORIENTED path on the double cover, so its coefficient is the negative - not
        # the same value, and not an unrelated one.
        other = coefficients(weber_voros_coefficient(p, m; order = 5, which = 1))
        for k in 1:2:5
            @test Complex{Float64}(other[k]) ≈ -Complex{Float64}(pred[k]) rtol = 1e-6
        end
    end

    @testset "ledger item 3: independent of the loop radius" begin
        # S_m dz has only even powers of the local coordinate t (z − z₀ = t²) on the
        # double cover, so the small loop carries no residue and δ cannot matter.
        p = weber_problem(1 // 1)
        m = weber_model(p)
        vals = [coefficients(weber_voros_coefficient(p, m; order = 3, loop_radius = δ))[1]
                for δ in (0.5, 1.0, 1.5)]
        @test Complex{Float64}(vals[1]) ≈ Complex{Float64}(vals[2]) rtol = 1e-12
        @test Complex{Float64}(vals[2]) ≈ Complex{Float64}(vals[3]) rtol = 1e-12
    end

    @testset "log Γ by resurgence: the functional equation" begin
        # No external Gamma anywhere - the value comes from Borel-summing the Weber
        # Voros coefficient. Γ(z+1) = zΓ(z) at z = ½ + ν is the sharpest test that
        # does not need one either.
        for ν in (3.0, 6.0, 12.0)
            lg = weber_log_gamma(ν; order = 11)
            lg1 = weber_log_gamma(ν + 1; order = 11)
            @test real(lg1 - lg) ≈ log(ν + 0.5) atol = 1e-10
            @test imag(lg) ≈ 0 atol = 1e-12
        end
    end

    @testset "log Γ by resurgence: the absolute anchor Γ(½) = √π" begin
        # The functional equation fixes the value only up to a constant; walking the
        # recurrence down from ν = 10 to ν = 0 pins it against √π.
        lg = weber_log_gamma(10.0; order = 13)
        for k in 9:-1:0
            lg -= log(k + 0.5)
        end
        @test real(lg) ≈ log(π) / 2 atol = 1e-12
        # and the connection constant √(2π)/Γ(½+ν) at ν = 0 is then √2
        @test abs(weber_connection(10.0; order = 13)) ≈
              sqrt(2π) / exp(real(weber_log_gamma(10.0; order = 13))) rtol = 1e-12
    end

    @testset "the recurrence shift reaches ν = 0 and the imaginary axis" begin
        # ledger item 7: the Borel singularities at 2πi k sit ON the Laplace ray that
        # ν = iε asks for, so the shift is not an optimisation - it is the only route.
        @test real(weber_log_gamma(0.0; order = 13)) ≈ log(π) / 2 atol = 1e-12
        @test imag(weber_log_gamma(0.0; order = 13)) ≈ 0 atol = 1e-14
        # ... and it changes nothing where the direct sum already worked
        @test weber_voros_sum(12.0; order = 11, shift = 0) ≈
              weber_voros_sum(12.0; order = 11) rtol = 1e-12
        # functional equation, now at small ν where only the shift can evaluate it
        for ν in (0.0, 0.25, 0.5im, -1.5im)
            lg, lg1 = weber_log_gamma(ν; order = 11), weber_log_gamma(ν + 1; order = 11)
            @test lg1 - lg ≈ log(ν + 0.5) atol = 1e-11
        end
    end

    @testset "ledger item 6: the log-ν branch flip on the imaginary axis" begin
        # Re(Voros sum) at ν = iε is elementary - −½log(1 + e^{−2π|ε|}) - and it is the
        # SAME formula for both signs of ε only because the principal log ν flips at
        # ε = 0. That cancellation is what makes the uniform quantization condition of
        # `quantization.jl` real above the barrier top as well as below.
        for ε in (-3.0, -1.0, -0.2, 0.0, 0.2, 1.0, 3.0)
            G = weber_voros_sum(complex(0.0, ε); order = 13)
            @test real(G) ≈ -log(1 + exp(-2π * abs(ε))) / 2 atol = 1e-12
            @test imag(G) ≈ -imag(weber_voros_sum(complex(0.0, -ε); order = 13)) atol = 1e-12
        end
        @test real(weber_voros_sum(0.0im; order = 13)) ≈ -log(2) / 2 atol = 1e-12
    end

    @testset "the barrier amplitude by resurgence" begin
        # |T| = e^{−πε/2}/|√(2π)/Γ(½+iε)| = (1 + e^{2πε})^{−1/2}: the package's own
        # Borel-summed Γ reproduces the parabolic-cylinder transmission amplitude,
        # including exactly at the barrier top (ε = 0, |T| = 1/√2).
        for ε in (-3.0, -1.0, -0.2, 0.0, 0.2, 1.0, 3.0)
            @test weber_barrier_amplitude(ε; order = 13) ≈ 1 / sqrt(1 + exp(2π * ε)) rtol = 1e-12
        end
        @test weber_barrier_amplitude(0.0; order = 13) ≈ 1 / sqrt(2) rtol = 1e-12
        # deep limit: the one-instanton weight e^{−πε}; over-barrier limit: full
        # transmission
        @test weber_barrier_amplitude(4.0; order = 13) ≈ exp(-4π) rtol = 1e-6
        @test weber_barrier_amplitude(-6.0; order = 13) ≈ 1 rtol = 1e-14
    end

    @testset "Borel singularities of the Stirling series sit at 2πi k" begin
        # The Gamma function's own resurgence, on a series whose singularity lattice
        # is known exactly - and a workout for the Borel-plane layer.
        S = weber_voros_series(; order = 21)
        Sf = FormalSeries(ComplexF64[ComplexF64(c) for c in coefficients(S)], :x;
                          power_offset = 1 // 1)
        B = borel(Sf)
        ps = sort(poles(pade(B; reduce = true)); by = abs)
        @test length(ps) ≥ 4
        for z in ps[1:4]
            @test abs(real(z)) < 1e-6 * (1 + abs(z))          # purely imaginary
            @test abs(imag(z) / (2π) - round(imag(z) / (2π))) < 1e-4
        end
        @test sort(abs.(imag.(ps[1:4])) / (2π))[[1, 3]] ≈ [1.0, 2.0] rtol = 1e-4
    end

    @testset "harmonic exactness coexists with a divergent Voros coefficient" begin
        # The CLOSED merging cycle of the Weber problem has all odd quantum periods
        # zero (the WKB series truncates - the harmonic oracle of test_voros.jl),
        # while the OPEN path to infinity is a factorially divergent series. Both are
        # true at once; that is exactly why the Voros coefficient needs the open path.
        p = weber_problem(1 // 1)
        m = weber_model(p)
        w = wkb_expansion(p; order = 5)
        closed = voros_symbol(w, encircling_contour(m.turning_points...))
        @test all(abs.(Complex{Float64}.(coefficients(quantum_series(closed)))) .< 1e-6)
        open_c = coefficients(weber_voros_coefficient(p, m; order = 5))
        @test abs(Complex{Float64}(open_c[1])) > 1e-3
        # Gevrey-1 growth: |c_n| ~ (2n−2)!/(2π)^{2n}, so ratios grow without bound
        ex = coefficients(weber_voros_series(; order = 21))
        r = [abs(Float64(ex[k + 2]) / Float64(ex[k])) for k in 1:2:17]
        @test issorted(r[3:end])
        @test r[end] > 10 * r[1]
    end
end
