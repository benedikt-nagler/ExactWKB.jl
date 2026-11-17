using Resurgence: FormalSeries, coefficients, power_offset, borel, pade

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
        # ε = 0. This test concerns only the Weber gamma expression.
        for ε in (-3.0, -1.0, -0.2, 0.0, 0.2, 1.0, 3.0)
            G = weber_voros_sum(complex(0.0, ε); order = 13)
            @test real(G) ≈ -log(1 + exp(-2π * abs(ε))) / 2 atol = 1e-12
            @test imag(G) ≈ -imag(weber_voros_sum(complex(0.0, -ε); order = 13)) atol = 1e-12
        end
        @test real(weber_voros_sum(0.0im; order = 13)) ≈ -log(2) / 2 atol = 1e-12
    end


    @testset "Borel singularities of the Stirling series sit at 2πi k" begin
        # The Gamma function's own resurgence, on a series whose singularity lattice
        # is known exactly - and a workout for the Borel-plane layer.
        S = weber_voros_series(; order = 21)
        Sf = FormalSeries(ComplexF64[ComplexF64(c) for c in coefficients(S)], :x;
                          power_offset = 1 // 1)
        B = borel(Sf)
        ps = sort(Resurgence.poles(pade(B; reduce = true)); by = abs)
        @test length(ps) ≥ 4
        for z in ps[1:4]
            @test abs(real(z)) < 1e-6 * (1 + abs(z))          # purely imaginary
            @test abs(imag(z) / (2π) - round(imag(z) / (2π))) < 1e-4
        end
        @test sort(abs.(imag.(ps[1:4])) / (2π))[[1, 3]] ≈ [1.0, 2.0] rtol = 1e-4
    end


end

@testset "Weber gamma against a separate Euler integral" begin
    import QuadGK
    for nu in (0.3,1.0,4.0,12.0)
        reference,_ = QuadGK.quadgk(t -> t^(nu-0.5)*exp(-t),0.0,Inf;rtol=1e-12)
        @test exp(weber_log_gamma(nu;order=13)) ≈ reference rtol=1e-10
    end
    for epsilon in (-2.0,-0.3,0.0,0.3,2.0)
        @test abs2(exp(weber_log_gamma(im*epsilon;order=13))) ≈ pi/cosh(pi*epsilon) rtol=1e-11
    end
end
