using ExactWKB
using ExactWKB: pair, theta
using QuadGK: quadgk
using Test

@testset "saddles" begin
    @testset "Airy has no saddles" begin
        airy = SchrodingerProblem([0, 1])          # one turning point → no pairs
        @test isempty(saddle_candidates(airy))
        @test isempty(saddles(airy))
    end

    @testset "double well: 3 BPS states" begin
        dw = SchrodingerProblem([3//4, 0, -2, 0, 1])   # (z²−1)² − 1/4
        cands = saddle_candidates(dw)
        @test length(cands) == 6                        # C(4,2) pairs

        # verification keeps only the three primitive saddles:
        #   inner pair (2,3) = instanton, outer-inner (1,2)/(3,4) = well cycles.
        sad = saddles(dw)
        @test length(sad) == 3
        @test Set(pair.(sad)) == Set([(1, 2), (2, 3), (3, 4)])

        # --- instanton (inner pair): real central charge, θ_c = 0 ---
        inst = only(filter(s -> pair(s) == (2, 3), sad))
        @test isapprox(theta(inst), 0.0; atol = 1e-6)
        @test abs(imag(central_charge(inst))) < 1e-6    # Z is real

        # independent oracle: |Z| = 2∫√Q over the real inner segment where Q ≥ 0
        a = sqrt(0.5)
        Qreal(z) = (z^2 - 1)^2 - 0.25
        Zref, _ = quadgk(z -> sqrt(Qreal(z)), -a, a; rtol = 1e-12)
        @test isapprox(mass(inst), 2 * abs(Zref); atol = 1e-6)

        # --- perturbative well cycles: purely imaginary Z, θ_c = π/2, equal mass ---
        cyc = filter(s -> pair(s) in ((1, 2), (3, 4)), sad)
        @test length(cyc) == 2
        for s in cyc
            @test isapprox(theta(s), π / 2; atol = 1e-6)
            @test abs(real(central_charge(s))) < 1e-6
        end
        @test isapprox(mass(cyc[1]), mass(cyc[2]); rtol = 1e-8)   # left/right symmetry
        @test mass(cyc[1]) < mass(inst)                          # cycle lighter than instanton
    end

    @testset "is_saddle agrees with the traced graph" begin
        dw = SchrodingerProblem([3//4, 0, -2, 0, 1])
        g0 = stokes_graph(dw; theta = 0.0)
        @test is_saddle(g0, (2, 3))          # instanton edge present at θ = 0
        @test is_saddle(g0, (3, 2))          # order-insensitive
        @test !is_saddle(g0, (1, 2))         # well cycle absent at θ = 0
    end

    @testset "topology jumps across the instanton's critical phase" begin
        dw = SchrodingerProblem([3//4, 0, -2, 0, 1])
        fam = stokes_graph_family(dw, [-0.05, 0.0, 0.05])
        @test isempty(edges(fam[1]))         # below θ_c: no finite edge
        @test edges(fam[2]) == [(2, 3)]      # at θ_c = 0: the instanton appears
        @test isempty(edges(fam[3]))         # above θ_c: gone again
        # the number of infinite rays jumps by 2 as the saddle forms
        @test n_infinite_lines(fam[2]) == n_infinite_lines(fam[1]) - 2
    end

    @testset "verify = false returns raw candidates" begin
        dw = SchrodingerProblem([3//4, 0, -2, 0, 1])
        @test length(saddles(dw; verify = false)) == 6
    end
end
