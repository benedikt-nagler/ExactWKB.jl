@testset "turning_points" begin
    @testset "Airy Q = z: one simple TP at 0" begin
        tps = turning_points(SchrodingerProblem([0, 1]))
        @test length(tps) == 1
        @test is_simple(tps[1])
        @test order(tps[1]) == 1
        @test abs(location(tps[1])) < 1e-30
    end

    @testset "harmonic ±√E to BigFloat precision" begin
        setprecision(256) do
            E = 2
            tps = turning_points(SchrodingerProblem([0, 0, 1]; energy = E))
            @test length(tps) == 2
            @test all(is_simple, tps)
            s = sqrt(BigFloat(E))
            # sorted by (real, imag): −√2 then +√2
            @test abs(location(tps[1]) - (-s)) < 1e-70
            @test abs(location(tps[2]) - s) < 1e-70
        end
    end

    @testset "double well → four simple TPs" begin
        # (z²−1)² − 1/4 = z⁴ − 2z² + 3/4, roots z = ±√(3/2), ±√(1/2)
        dwell = SchrodingerProblem([3//4, 0, -2, 0, 1])
        tps = turning_points(dwell)
        @test length(tps) == 4
        @test all(is_simple, tps)
        locs = sort(real.(location.(tps)))
        @test locs ≈ [-sqrt(3/2), -sqrt(1/2), sqrt(1/2), sqrt(3/2)]
        @test all(x -> abs(imag(x)) < 1e-30, location.(tps))
    end

    @testset "multiplicity: Q = z² and (z−1)²(z+2)" begin
        double = turning_points(SchrodingerProblem([0, 0, 1]))  # z²
        @test length(double) == 1
        @test order(double[1]) == 2
        @test !is_simple(double[1])
        @test abs(location(double[1])) < 1e-20

        # (z−1)²(z+2) = z³ − 3z + 2
        p = SchrodingerProblem([2, -3, 0, 1])
        tps = turning_points(p)
        @test length(tps) == 2
        byloc = sort(tps; by = t -> real(location(t)))
        @test real(location(byloc[1])) ≈ -2 && order(byloc[1]) == 1   # z = −2 simple
        @test real(location(byloc[2])) ≈ 1  && order(byloc[2]) == 2   # z = 1 double
        @test simple_turning_points(p) == [byloc[1]]
    end

    @testset "Float64 potential stays Float64" begin
        tps = turning_points(SchrodingerProblem([0.0, 0.0, 1.0]; energy = 1.0))
        @test eltype(location.(tps)) == Complex{Float64}
        @test sort(real.(location.(tps))) ≈ [-1.0, 1.0]
    end
end
