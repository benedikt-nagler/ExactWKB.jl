# M5 flagship demonstration: spectra from exact WKB — energy-dependent quantum
# periods, exact (Voros) quantization conditions, median-summed eigenvalues, and the
# double-well ZJJ/Dunne–Ünsal resurgent structure.
#
# Three acts:
#   1. Quartic oscillator V = z² + z⁴: `wkb_eigenvalue` vs dense diagonalization —
#      a digit-count table of the median-summed exact-WKB spectrum.
#   2. Symmetric double well V = (z²−1)²: the even/odd splitting vs ħ on a log plot
#      against the one-instanton line  ~ e^{−A_cl/2ħ}.
#   3. The Dunne–Ünsal P/NP relation verified order by order (c = 3/2 for this
#      potential), and the Borel plane of the barrier Voros symbol — the A/−A
#      resonant ray that Resurgence M6b's log sectors were built for.
#
# Run:  julia --project=test examples/m5_flagship.jl

using ExactWKB
import Resurgence
using CairoMakie
using Printf

# harmonic-basis diagonalization oracle (test-side helper, LinearAlgebra stdlib only)
include(joinpath(@__DIR__, "..", "test", "diagonalization.jl"))

# ── 1. Quartic spectrum: exact WKB vs diagonalization ───────────────────────────────
println("── Quartic oscillator V = z² + z⁴ ", "─"^42)
quart = SchrodingerProblem([0, 0, 1, 0, 1])
ħ = 0.05
E_diag = setprecision(192) do
    diagonalization_eigenvalues(quart, ħ; nev = 4, N = 120)
end
println("ħ = $ħ; exact quantization condition: Im log V_well + 2π(n+½) = 0\n")
println("  n   E (exact WKB)             E (diagonalization)       digits")
for n in 0:3
    E = setprecision(192) do
        wkb_eigenvalue(quart, n, big(ħ); order = 16, rtol = 1e-15, quad_rtol = 1e-25)
    end
    digits = -log10(Float64(abs(E - E_diag[n + 1]) / abs(E_diag[n + 1])))
    @printf("  %d   %.18f     %.18f     %.1f\n", n, Float64(E),
            Float64(E_diag[n + 1]), digits)
end

# ── 2. Double-well splitting vs the one-instanton line ─────────────────────────────
println("\n── Double well V = (z²−1)²: even/odd splitting ", "─"^30)
dwell = SchrodingerProblem([1, 0, -2, 0, 1])
hs = [0.08, 0.10, 0.12, 0.15, 0.18]
splittings = Float64[]
instanton = Float64[]
for h in hs
    ΔE, pred = setprecision(160) do
        Δ = energy_splitting(dwell, 0, big(h); order = 14, rtol = 1e-12,
                             quad_rtol = 1e-20)
        # one-instanton prediction from the periods themselves: the parity condition
        # cos φ = ±√(V_A/(1+V_A)) linearizes to ΔE ≈ 2√V_A / |dφ/dE|, φ = πB
        E0 = wkb_eigenvalue(dwell, 0, big(h); parity = +1, order = 14, rtol = 1e-10,
                            quad_rtol = 1e-18)
        δ = big(1) / 1000
        A = instanton_a(dwell, E0; order = 8, quad_rtol = 1e-18)
        Bp = perturbative_b(dwell, E0 + δ; order = 8, quad_rtol = 1e-18)
        Bm = perturbative_b(dwell, E0 - δ; order = 8, quad_rtol = 1e-18)
        dφdE = big(π) * real(Resurgence.evaluate(Bp - Bm, big(h))) / (2δ)
        VA = exp(-real(Resurgence.evaluate(A, big(h))))
        # ΔE ≈ Δφ/(dφ/dE), Δφ = 2·(½√V_A) = √V_A  (the ½ is the parity condition's)
        (Float64(Δ), Float64(sqrt(VA) / abs(dφdE)))
    end
    push!(splittings, ΔE)
    push!(instanton, pred)
    @printf("  ħ = %.2f:  ΔE = %.6e   one-instanton ≈ %.6e\n", h, ΔE, pred)
end
fig = Figure(size = (640, 420))
ax = Axis(fig[1, 1]; xlabel = "1/ħ", ylabel = "ΔE₀", yscale = log10,
          title = "Double-well splitting: exact WKB vs one-instanton asymptotics")
scatter!(ax, 1 ./ hs, splittings; label = "energy_splitting (median-summed)")
lines!(ax, 1 ./ hs, instanton; color = :orangered,
       label = "√V_A / (dφ/dE)  (one-instanton)")
axislegend(ax; position = :rt)
out1 = joinpath(@__DIR__, "m5_splitting.png")
CairoMakie.save(out1, fig)
println("Wrote ", out1)

# ── 3. Dunne–Ünsal P/NP relation + the resonant Borel plane ─────────────────────────
println("\n── Dunne–Ünsal relation 1/(∂B/∂E) = −c·ħ³·∂A/∂ħ|_B ", "─"^25)
r = setprecision(160) do
    verify_zjj(dwell, big(1) / 2; order = 9, quad_rtol = 1e-20)
end
@printf("fitted c = %.10f   (exact: 3/2 for V = (z²−1)²)\n", Float64(real(r.c)))
for (p, res) in zip(r.orders, r.residuals)
    @printf("  ħ^%-2d residual: %.2e\n", p, Float64(res))
end

# the barrier Voros symbol's Borel plane: singularities on the A/−A resonant ray
# (the odd-only Voros series needs the reduced Padé, then plot the approximant)
vsA = quantum_period(dwell, :barrier, 0.5; order = 12)
rA = Resurgence.pade(Resurgence.borel(quantum_series(vsA)); reduce = true)
figB = Resurgence.plot_borel_plane(rA; rays = [0.0])
out2 = joinpath(@__DIR__, "m5_borel_barrier.png")
CairoMakie.save(out2, figB)
println("\nBarrier symbol Borel plane (A/−A poles on the real axis): wrote ", out2)
println("\nM5 flagship complete.")
