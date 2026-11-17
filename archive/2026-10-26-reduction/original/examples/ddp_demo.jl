# DDP flagship demonstration: the cubic potential Q = z³ − z (type A₂) run through the
# convention-pinned DDP layer - saddles and their critical phases, Voros symbols with
# Borel poles on the Z-lattice, the measured lateral-sum jump vs the DDP prediction
# (Stokes constant = an exact integer = the intersection pairing), the jump as a
# cluster y-mutation, and the pentagon / DT closure of the chamber.
#
# Run:  julia --project=test examples/ddp_demo.jl

using ExactWKB
import Resurgence
import ClusterAlgebras
using CairoMakie

prob = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])       # Q = z³ − z
println(prob)
tps = simple_turning_points(prob)
println("\nTurning points: ", [location(t) for t in tps])

# ── 1. BPS saddles and their walls ─────────────────────────────────────────────────
cands = saddle_candidates(prob)
sads = saddles(prob)
println("\nSaddle candidates (Z = 2∫√Q, θ_c = arg Z mod π):")
for s in cands
    confirmed = any(t -> t.pair == s.pair, sads) ? "confirmed" : "not a saddle here"
    println("  ", s.pair, "  Z = ", round(central_charge(s); sigdigits = 6),
            "  θ_c = ", round(ExactWKB.theta(s); sigdigits = 6), "  (", confirmed, ")")
end
Z = Dict(s.pair => central_charge(s) for s in cands)
println("Charge additivity: Z₁₃ − (Z₁₂ − Z₂₃) = ",
        round(abs(Z[(1, 3)] - (Z[(1, 2)] - Z[(2, 3)])); sigdigits = 3),
        "  → two-state chamber (the composite is not a saddle at its phase)")

# ── 2. Voros symbols, decaying orientation, Borel poles on the Z-lattice ───────────
w = wkb_expansion(prob; order = 12)
c1 = reverse(encircling_contour(tps[1], tps[2]; margin = 0.4))   # Z₁ < 0: decaying
c2 = encircling_contour(tps[2], tps[3]; margin = 0.4)            # Z₂ ∈ i𝐑₋: decaying
v1 = voros_symbol(w, c1)
v2 = voros_symbol(w, c2)
println("\nVoros symbols (decaying orientation):")
println("  Z(γ₁) = ", round(classical_period(v1); sigdigits = 6),
        "   wall θ_c = 0")
println("  Z(γ₂) = ", round(classical_period(v2); sigdigits = 6),
        "   wall θ_c = π/2")
B2 = Resurgence.borel(quantum_series(v2))
r2 = Resurgence.pade(B2; reduce = true)
println("  Borel poles of γ₂'s series: ",
        round.(ComplexF64.(Resurgence.poles(r2)); sigdigits = 4))
println("  → the pole on the θ = 0 ray sits at ζ = −Z(γ₁) = ",
        round(-classical_period(v1); sigdigits = 6))

# ── 3. Figure 1: the saddle appearing and dying across the wall ────────────────────
thetas = [-0.15, 0.0, 0.15]
family = stokes_graph_family(prob, thetas)
fig1 = plot_stokes_graph(family)
out1 = joinpath(@__DIR__, "ddp_wall_crossing.png")
CairoMakie.save(out1, fig1)
println("\nWrote ", out1, "  (Stokes graphs at θ = ", thetas,
        ": the finite γ₁ edge exists only on the wall)")

# ── 4. The DDP jump: Stokes constant = exact integer ───────────────────────────────
κ = intersection_pairing(prob, c2, c1)
println("\nIntersection pairing ⟨γ₂, γ₁⟩ = ", κ)
Zwall = abs(classical_period(v1))
println("DDP jump of log V(γ₂) across θ = 0 vs κ·log(1 + V(γ₁)):")
xs = [4.0, 5.0, 6.0, 7.0, 8.0]
jumps = Float64[]
preds = Float64[]
for x in xs
    ħ = Zwall / x
    r = verify_ddp(v2, v1, κ, ħ; theta = 0.0)
    push!(jumps, abs(r.jump_measured))
    push!(preds, abs(r.jump_predicted))
    println("  |Z|/ħ = ", x, ":  measured = ", round(r.jump_measured; sigdigits = 5),
            "  predicted = ", round(r.jump_predicted; sigdigits = 5),
            "  κ_measured = ", round(real(r.kappa_measured); sigdigits = 4),
            "  → ", round(Int, real(r.kappa_measured)))
end

# ── 5. The jump as a cluster y-mutation ─────────────────────────────────────────────
P = [0 -κ; κ 0]                                        # P[i,j] = ⟨γ_i, γ_j⟩
seed = ddp_seed(P)
println("\nCluster seed from the pairing (B = −P): type ",
        ClusterAlgebras.cartan_type(ClusterAlgebras.Quiver(-P)))
ħ = Zwall / 5
res = verify_ddp_mutation([v1, v2], seed, 1, ħ; theta = 0.0)
println("Crossing the γ₁ wall = y-mutation μ₁: max residual = ",
        round(res.max_residual; sigdigits = 3))
ms = ClusterAlgebras.mutate(seed, 1)
println("Tropical shadow (c-vector mutation): ",
        ClusterAlgebras.cvectors(seed), " → ", ClusterAlgebras.cvectors(ms))

# ── 6. Pentagon / DT closure ────────────────────────────────────────────────────────
mgs = ClusterAlgebras.maximal_green_sequences(seed)
println("\nMaximal green sequences of the seed: ", mgs)
chamber = [[0, 1], [1, 0]]                             # θ-decreasing saddle charges
seq = mgs[findfirst(s -> ClusterAlgebras.ordered_c_vectors(seed, s) == chamber, mgs)]
println("Chamber (θ-decreasing charges ", chamber, ") = MGS ", seq)
word = ClusterAlgebras.quantum_dilog_word(seed, seq)
println("DT invariants Ω(γ): ", ClusterAlgebras.omega(word), "  (all 1)")
println("Zamolodchikov / pentagon periodicity: y-system period = ",
        ClusterAlgebras.y_system(seed).period, "  (= h + 2 = 5)")
dt = ClusterAlgebras.dt_transformation(seed)
println("DT transformation: sequence ", dt.sequence, ", σ = ", dt.sigma,
        ", tropical y_images = y ↦ y⁻¹")

# ── 7. Figure 2: measured jump vs prediction over ħ ─────────────────────────────────
fig2 = Figure(size = (640, 420))
ax = Axis(fig2[1, 1]; xlabel = "|Z(γ₁)| / ħ", ylabel = "|jump of log V(γ₂)|",
          yscale = log10, title = "DDP jump across θ = 0 (cubic, A₂)")
scatter!(ax, xs, jumps; label = "measured  s⁺ − s⁻", markersize = 12)
lines!(ax, xs, preds; label = "predicted  κ·log(1 + V(γ₁))", color = :orangered)
axislegend(ax; position = :rt)
out2 = joinpath(@__DIR__, "ddp_jump_vs_hbar.png")
CairoMakie.save(out2, fig2)
println("\nWrote ", out2)
