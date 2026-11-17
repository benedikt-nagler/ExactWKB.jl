@testset "potentials" begin
    @testset "q_taylor_at: exact Taylor shift" begin
        # Q = (z−1)²(z+2) = z³ − 3z + 2, expanded about z₀
        p = SchrodingerProblem([2, -3, 0, 1])
        @test q_taylor_at(p, 0) == [2, -3, 0, 1]                  # already at 0
        @test q_taylor_at(p, 1) == [0, 0, 3, 1]                   # double zero at 1
        @test q_taylor_at(p, -2) == [0, 9, -6, 1]                 # simple zero at −2
        @test eltype(q_taylor_at(p, 1 // 2)) <: Rational          # stays exact
        # t[1] = Q(z₀), t[2] = Q′(z₀) reproduce the existing accessors
        for z in (0.3, -1.7, 2.0)
            t = q_taylor_at(p, z)
            @test t[1] ≈ p(z)
            @test t[2] ≈ q_derivative_at(p, z)
        end
        # length is degree + 1, and the top coefficient is z₀-independent
        @test length(q_taylor_at(p, 5)) == degree(p) + 1
        @test q_taylor_at(p, 5)[end] == last(q_coefficients(p))
    end

    @testset "construction & Q = V − E" begin
        airy = SchrodingerProblem([0, 1])                 # Q = z
        @test q_coefficients(airy) == [0, 1]
        @test degree(airy) == 1
        @test variable(airy) === :z
        @test energy(airy) == 0

        harm = SchrodingerProblem([0, 0, 1]; energy = 1)  # Q = z² − 1
        @test q_coefficients(harm) == [-1, 0, 1]
        @test degree(harm) == 2
        @test energy(harm) == 1

        # exact rational coefficients stay exact
        dwell = SchrodingerProblem([3//4, 0, -2, 0, 1])   # (z²−1)² − 1/4
        @test eltype(q_coefficients(dwell)) == Rational{Int}
        @test degree(dwell) == 4
    end

    @testset "invalid potentials" begin
        @test_throws InvalidPotential SchrodingerProblem(Int[])
        # a genuinely constant Q (all non-constant coefficients vanish)
        @test_throws InvalidPotential SchrodingerProblem([5])
        @test_throws InvalidPotential SchrodingerProblem([5, 0, 0]; energy = 5)
        # Q = V − E constant even though V is not (only the constant term differs)
        @test SchrodingerProblem([1, 2]; energy = 2) isa SchrodingerProblem  # Q = −1 + 2z ok
    end

    @testset "evaluation & derivative" begin
        harm = SchrodingerProblem([0, 0, 1]; energy = 1)  # Q = z² − 1, Q′ = 2z
        @test harm(2.0) == 3.0
        @test harm(0) == -1
        @test q_derivative_at(harm, 3.0) == 6.0
        @test q_derivative_at(harm, 0) == 0

        # promotion: exact coeffs, complex argument
        v = harm(1.0 + 1.0im)                              # (1+i)² − 1 = -1 + 2i
        @test v ≈ -1 + 2im
        @test v isa Complex{Float64}

        # exact evaluation stays exact
        @test harm(2 // 1) == 3 // 1
        @test harm(2 // 1) isa Rational
    end

    @testset "with_energy" begin
        harm = SchrodingerProblem([0, 0, 1])
        h1 = with_energy(harm, 4)
        @test energy(h1) == 4
        @test energy(harm) == 0                            # original unchanged
        @test q_coefficients(h1) == [-4, 0, 1]
        @test h1 == SchrodingerProblem([0, 0, 1]; energy = 4)
    end
end
