# M4 flagship demonstration: the full cluster bridge, no hand-built input.
#
# Given a polynomial potential, the bridge computes everything the DDP layer used to
# be fed by hand: Stokes graph → ideal triangulation → quiver (`ideal_triangulation`,
# `triangulation_quiver`), the cycle basis with central charges and intersection
# pairing (`charge_basis`), the cluster seed (`bridge_seed`), and the BPS spectrum
# from the greedy phase-ordered maximal green sequence (`bps_spectrum`). The flagship
# verification: Borel-plane singularities of the numerically summed Voros symbols
# land on the exact central-charge lattice Z_γ, and the Stokes constants are the
# exact integers Ω(γ)·⟨γ,γ'⟩ — cluster data.
#
# Run:  julia --project=test examples/m4_flagship.jl

using ExactWKB
import Resurgence
import ClusterAlgebras
using CairoMakie

# ── 1. Cubic (A₂): graph → triangulation → quiver ───────────────────────────────────
cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])       # Q = z³ − z
println(cubic)
g = stokes_graph(cubic; theta = 0.3)
t = ideal_triangulation(g)
println("\n", g)
show(stdout, MIME"text/plain"(), t); println()
q = triangulation_quiver(t)
println("Triangulation quiver: B = ", q.B, ",  type ",
        ClusterAlgebras.cartan_type(q), "  (pentagon → A₂)")
fig1 = plot_triangulation(g, t)
out1 = joinpath(@__DIR__, "m4_cubic_triangulation.png")
CairoMakie.save(out1, fig1)
println("Wrote ", out1)

# ── 2. Automatic charge basis (what ddp_demo.jl built by hand) ──────────────────────
cb = charge_basis(cubic, t; margin = 0.4)
show(stdout, MIME"text/plain"(), cb); println()
println("Keystone: |P| == |B| : ", abs.(cb.pairing) == abs.(q.B),
        "   (hand-built ddp_demo values: Z₁ ≈ -0.958512, Z₂ ≈ -0.958512i, P = [0 -1; 1 0])")
seed = bridge_seed(cb)
println("bridge_seed: B = -P, type ",
        ClusterAlgebras.cartan_type(ClusterAlgebras.Quiver(-cb.pairing)))

# ── 3. The BPS spectrum from the greedy maximal green sequence ──────────────────────
sp = bps_spectrum(cubic; theta = 0.3, margin = 0.4)
show(stdout, MIME"text/plain"(), sp); println()
println("Physical MGS: ", sp.sequence, "  (charges in θ-decreasing wall order,")
println("verified against the traced saddles — ledger item 6)")
println("Pentagon periodicity: y-system period = ",
        ClusterAlgebras.y_system(seed).period, "  (= h + 2 = 5)")

# ── 4. Flagship verification: Borel poles on the Z-lattice, integer Stokes constants ─
w = wkb_expansion(cubic; order = 12)
vs = [voros_symbol(w, c) for c in sp.basis.contours]
k1 = findfirst(Z -> abs(imag(Z)) < 1e-8, central_charges(sp.basis))  # the θ_c = 0 wall
k2 = 3 - k1
println("\nBorel poles of the Voros symbols vs the BPS lattice:")
lattice = vcat([central_charge(s) for s in sp.states],
               [-central_charge(s) for s in sp.states])
polesets = [ComplexF64.(Resurgence.poles(Resurgence.pade(
                Resurgence.borel(quantum_series(v)); reduce = true))) for v in vs]
for (j, ps) in enumerate(polesets)
    println("  γ", sp.basis.triangulation.diagonal_tp_pair[j], " poles: ",
            round.(ps; sigdigits = 4))
end
println("  BPS lattice ±Z_γ: ", round.(ComplexF64.(lattice); sigdigits = 4))
println("  → γ", sp.basis.triangulation.diagonal_tp_pair[k2],
        "'s pole on the θ = 0 ray sits at ζ = -Z(γ",
        sp.basis.triangulation.diagonal_tp_pair[k1], ") = ",
        round(ComplexF64(-central_charges(sp.basis)[k1]); sigdigits = 6))

κ = cb.pairing[k2, k1]
Zwall = abs(central_charges(sp.basis)[k1])
println("\nStokes constants across θ = 0 (exact integer = Ω(γ')·⟨γ,γ'⟩ = ", κ, "):")
xs = [4.0, 5.0, 6.5, 8.0]
for x in xs
    r = verify_ddp(vs[k2], vs[k1], κ, Zwall / x; theta = 0.0)
    println("  |Z|/ħ = ", x, ":  κ_measured = ",
            round(real(r.kappa_measured); sigdigits = 4), "  → ",
            round(Int, real(r.kappa_measured)),
            "   (relative residual ", round(r.relative_residual; sigdigits = 2), ")")
end
res = verify_ddp_mutation(vs, seed, k1, Zwall / 5; theta = 0.0)
println("Wall-crossing = y-mutation μ_", k1, ": max residual = ",
        round(res.max_residual; sigdigits = 3), "  (fully automatic input)")

fig2 = Figure(size = (640, 480))
ax2 = Axis(fig2[1, 1]; xlabel = "Re ζ", ylabel = "Im ζ", aspect = DataAspect(),
           title = "Borel plane: measured poles on the BPS lattice (cubic, A₂)")
for (j, ps) in enumerate(polesets)
    kept = filter(z -> abs(z) < 3.5, ps)
    scatter!(ax2, real.(kept), imag.(kept); markersize = 10,
             label = "poles of V(γ$(sp.basis.triangulation.diagonal_tp_pair[j]))")
end
scatter!(ax2, real.(ComplexF64.(lattice)), imag.(ComplexF64.(lattice));
         marker = :xcross, markersize = 16, color = :black, label = "±Z_γ (exact)")
axislegend(ax2; position = :rt, framevisible = false)
out2 = joinpath(@__DIR__, "m4_borel_lattice.png")
CairoMakie.save(out2, fig2)
println("Wrote ", out2)

# ── 5. Quartic double well (A₃) ─────────────────────────────────────────────────────
quartic = SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0])
println("\n", quartic)
spq = bps_spectrum(quartic)
show(stdout, MIME"text/plain"(), spq); println()
gq = stokes_graph(quartic; theta = spq.theta0)
tq = ideal_triangulation(gq)
println("Hexagon triangulation, quiver type ",
        ClusterAlgebras.cartan_type(triangulation_quiver(tq)))
seedq = bridge_seed(spq.basis)
println("Zamolodchikov periodicity: ", ClusterAlgebras.y_system(seedq).period,
        "  (= h + 2 = 6)")
dtq = ClusterAlgebras.dt_transformation(seedq)
Cf = ClusterAlgebras.cmatrix(dtq.seed)
println("DT closure: C_final[σ(j), j] = -1 for all j : ",
        all(Cf[dtq.sigma[j], j] == -1 for j in 1:3))
fig3 = plot_triangulation(gq, tq)
out3 = joinpath(@__DIR__, "m4_quartic_triangulation.png")
CairoMakie.save(out3, fig3)
println("Wrote ", out3)

# ── 6. Quintic (A₄) — the stress test beyond the hand-worked examples ───────────────
quintic = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])
println("\n", quintic)
spx = bps_spectrum(quintic)                      # finds the physical (top) chamber
show(stdout, MIME"text/plain"(), spx); println()
println("Chamber θ₀ = ", round(spx.theta0; sigdigits = 4),
        ",  states = ", n_states(spx), " = confirmed saddles = ",
        length(saddles(quintic)))
gx = stokes_graph(quintic; theta = spx.theta0)
tx = ideal_triangulation(gx)
println("Heptagon triangulation, quiver type ",
        ClusterAlgebras.cartan_type(triangulation_quiver(tx)))
fig4 = plot_triangulation(gx, tx)
out4 = joinpath(@__DIR__, "m4_quintic_triangulation.png")
CairoMakie.save(out4, fig4)
println("Wrote ", out4)
println("(M4 finding: chambers whose triangulation has a fully internal triangle, or")
println(" whose sweep fails the saddle gate, need the deferred signed-flip layer —")
println(" bps_spectrum searches the wall-free gaps for a chamber that passes.)")

# ── 7. Wall-crossing finale: flip = mutation, spectrum chamber-independent ──────────
thetas = [-0.15, 0.0, 0.15]
family = stokes_graph_family(cubic, thetas)
fig5 = plot_stokes_graph(family)
out5 = joinpath(@__DIR__, "m4_wall_crossing.png")
CairoMakie.save(out5, fig5)
println("\nWrote ", out5)
tp = ideal_triangulation(stokes_graph(cubic; theta = 0.15))
tm = ideal_triangulation(stokes_graph(cubic; theta = -0.15))
k = findfirst(==((1, 2)), tp.diagonal_tp_pair)
println("Crossing the (1,2) wall flips its diagonal: μ_", k, "(B₊) == B₋ : ",
        ClusterAlgebras.mutate(triangulation_quiver(tp), k).B ==
        triangulation_quiver(tm).B)
spa = bps_spectrum(cubic; theta = 0.3)
spb = bps_spectrum(cubic; theta = 2.0)
println("Spectrum from two different chambers agrees (masses): ",
        sort([mass(s) for s in spa.states]) ≈ sort([mass(s) for s in spb.states]))
println("\nM4 bridge complete: Stokes graph → triangulation → quiver → BPS spectrum,")
println("Borel singularities = Z_γ, Stokes constants = integer DT data. No hand input.")
