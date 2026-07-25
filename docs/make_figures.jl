# Regenerates the figures used by the documentation into docs/src/assets/.
#
# Not part of the Documenter build (which stays free of a plotting backend) - run it by
# hand when a figure's underlying numbers change:
#
#   julia --project=test docs/make_figures.jl
#
# The test environment is used because it already carries CairoMakie and the three
# packages by path.

using ExactWKB
import Resurgence
import ClusterAlgebras
using CairoMakie

const ASSETS = joinpath(@__DIR__, "src", "assets")
mkpath(ASSETS)
save_fig(name, fig) = (CairoMakie.save(joinpath(ASSETS, name), fig);
                       println("wrote ", name))

# ─────────────────────────────────────────────────────────────────────────────────────
# Tutorial 1: the cubic
# ─────────────────────────────────────────────────────────────────────────────────────

cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
g = stokes_graph(cubic; theta = 0.3)
t = ideal_triangulation(g)
cb = charge_basis(cubic, t; margin = 0.4)
sp = bps_spectrum(cubic; theta = 0.3, margin = 0.4)

save_fig("cubic_stokes.png", plot_stokes_graph(g))
save_fig("cubic_triangulation.png", plot_triangulation(g, t))
save_fig("cubic_wall_crossing.png",
         plot_stokes_graph(stokes_graph_family(cubic, [-0.15, 0.0, 0.15])))

# Borel poles of the two Voros symbols against the exact central-charge lattice.
# Exact (rational) input at 192-bit precision and WKB order 20 - the leading Padé pole
# then sits on |Z| to 3e-4. Only the two leading poles of each symbol are drawn; the
# rest are the approximant's rendering of the further-out Borel structure.
let exact_cubic = SchrodingerProblem([0 // 1, -1 // 1, 0 // 1, 1 // 1])
    leads = setprecision(192) do
        etps = turning_points(exact_cubic)
        w20 = wkb_expansion(exact_cubic; order = 20)
        map(((1, 2), (2, 3))) do (i, j)
            v = voros_symbol(w20, encircling_contour(etps[i], etps[j]))
            ps = ComplexF64.(Resurgence.poles(Resurgence.pade(
                     Resurgence.borel(quantum_series(v)); reduce = true)))
            sort(ps; by = abs)[1:2]
        end
    end
    lattice = ComplexF64.(vcat([central_charge(s) for s in sp.states],
                               [-central_charge(s) for s in sp.states]))
    labels = ["leading poles of the γ(1,2) symbol", "leading poles of the γ(2,3) symbol"]

    fig = Figure(size = (960, 480))
    ax = Axis(fig[1, 1]; xlabel = "Re ζ", ylabel = "Im ζ", aspect = DataAspect(),
              title = "Borel plane: measured singularities on the BPS lattice")
    scatter!(ax, real.(lattice), imag.(lattice); marker = :xcross, markersize = 20,
             color = :black, label = "± Z_γ  (exact)")
    for (lead, lab) in zip(leads, labels)
        scatter!(ax, real.(lead), imag.(lead); markersize = 11, label = lab)
    end
    ylims!(ax, -1.5, 1.5)
    axislegend(ax; position = :lb, framevisible = false, labelsize = 11)

    axz = Axis(fig[1, 2]; xlabel = "Re ζ", ylabel = "Im ζ",
               xticks = [-4e-4, 0.0, 4e-4],
               title = "zoom on ζ = +i|Z|")
    scatter!(axz, [0.0], [abs(lattice[1])]; marker = :xcross, markersize = 22, color = :black)
    zoomed = filter(z -> abs(z - im * abs(lattice[1])) < 5e-3, vcat(leads...))
    scatter!(axz, real.(zoomed), imag.(zoomed); markersize = 13, color = Cycled(1))
    xlims!(axz, -6e-4, 6e-4); ylims!(axz, abs(lattice[1]) - 2e-4, abs(lattice[1]) + 5e-4)
    save_fig("cubic_borel_lattice.png", fig)
end

# The thermodynamic-Bethe-ansatz solution: log X along the rapidity line, with the
# infrared plateau at the golden ratio
let sol = solve_tba(sp)
    fig = Figure(size = (680, 420))
    ax = Axis(fig[1, 1]; xlabel = "rapidity θ", ylabel = "X_γ(θ)",
              title = "Thermodynamic Bethe ansatz: the cubic's two Y-functions")
    for a in 1:n_states(sp)
        lines!(ax, sol.grid, real.(exp.(sol.gplus[:, a]));
               linewidth = a == 1 ? 4 : 2, linestyle = a == 1 ? :solid : :dash,
               label = "X_γ$(a)  (the two coincide: equal |Z|)")
    end
    hlines!(ax, [(1 + sqrt(5)) / 2]; color = :black, linestyle = :dot,
            label = "golden ratio (constant Y-system)")
    axislegend(ax; position = :rt, framevisible = false)
    save_fig("cubic_tba.png", fig)
end

# ─────────────────────────────────────────────────────────────────────────────────────
# Tutorial 2: the double well
# ─────────────────────────────────────────────────────────────────────────────────────

dwell = SchrodingerProblem([1.0, 0.0, -2.0, 0.0, 1.0])

save_fig("dwell_stokes.png", plot_stokes_graph(stokes_graph(with_energy(dwell, 0.25);
                                                            theta = 0.0)))

# The potential with the ground doublet drawn in
let ħ = 0.1
    E_even = wkb_eigenvalue(dwell, 0, ħ; parity = +1, order = 8)
    E_odd = wkb_eigenvalue(dwell, 0, ħ; parity = -1, order = 8)
    zs = range(-1.6, 1.6; length = 400)
    fig = Figure(size = (680, 430))
    ax = Axis(fig[1, 1]; xlabel = "z", ylabel = "V(z)",
              title = "Double well V = (z²−1)²: the ground doublet at ħ = 0.1")
    lines!(ax, zs, [(z^2 - 1)^2 for z in zs]; color = :black, linewidth = 2)
    hlines!(ax, [E_even]; color = :seagreen, linewidth = 2, label = "even, E = $(round(E_even; sigdigits = 8))")
    hlines!(ax, [E_odd]; color = :crimson, linewidth = 2, linestyle = :dash,
            label = "odd,  E = $(round(E_odd; sigdigits = 8))")
    ylims!(ax, -0.05, 1.1)
    text!(ax, 0.0, 0.30; text = "ΔE = $(round(E_odd - E_even; sigdigits = 6))",
          align = (:center, :center), fontsize = 14)
    axislegend(ax; position = :ct, framevisible = false)
    save_fig("dwell_levels.png", fig)
end

# Splitting against the one-instanton line
let hs = [0.08, 0.10, 0.12, 0.15, 0.18]
    splittings = Float64[]
    instanton = Float64[]
    for h in hs
        push!(splittings, energy_splitting(dwell, 0, h; order = 8))
        E0 = wkb_eigenvalue(dwell, 0, h; parity = +1, order = 8)
        δ = 1e-4
        A = instanton_a(dwell, E0; order = 6)
        dφdE = π * real(Resurgence.evaluate(perturbative_b(dwell, E0 + δ; order = 6) -
                                            perturbative_b(dwell, E0 - δ; order = 6), h)) / (2δ)
        push!(instanton, sqrt(exp(-real(Resurgence.evaluate(A, h)))) / abs(dφdE))
    end
    fig = Figure(size = (660, 430))
    ax = Axis(fig[1, 1]; xlabel = "1/ħ", ylabel = "ΔE₀", yscale = log10,
              title = "Level splitting: exact WKB versus one-instanton asymptotics")
    scatter!(ax, 1 ./ hs, splittings; markersize = 12, label = "energy_splitting (resummed)")
    lines!(ax, 1 ./ hs, instanton; color = :orangered, label = "one-instanton estimate")
    axislegend(ax; position = :rt, framevisible = false)
    save_fig("dwell_splitting.png", fig)
end

# Borel plane of the barrier symbol: singularities on the well cycle's central charge
let cyc = spectral_cycles(dwell, 0.5)
    vw = quantum_period(dwell, cyc[3], 0.5; order = 8)
    vb = quantum_period(dwell, cyc[2], 0.5; order = 12)
    ps = ComplexF64.(Resurgence.poles(Resurgence.pade(
             Resurgence.borel(quantum_series(vb)); reduce = true)))
    Zw = ComplexF64(classical_period(vw))
    fig = Figure(size = (600, 500))
    ax = Axis(fig[1, 1]; xlabel = "Re ζ", ylabel = "Im ζ", aspect = DataAspect(),
              title = "Borel plane of the barrier symbol (E = 0.5)")
    scatter!(ax, real.(ps), imag.(ps); markersize = 12, label = "measured singularities")
    scatter!(ax, [real(Zw), -real(Zw)], [imag(Zw), -imag(Zw)]; marker = :xcross,
             markersize = 17, color = :black, label = "± Z of the well cycle (exact)")
    axislegend(ax; position = :rt, framevisible = false)
    save_fig("dwell_borel.png", fig)
end

# ─────────────────────────────────────────────────────────────────────────────────────
# Tutorial 3: Seiberg-Witten SU(2)
# ─────────────────────────────────────────────────────────────────────────────────────

sw = SeibergWittenSU2()

# The u-plane: wall of marginal stability, singular points, chambers, monodromy loops
let wall = ms_wall(sw; n = 96)
    closed = vcat(wall, wall[1:1])
    fig = Figure(size = (640, 560))
    ax = Axis(fig[1, 1]; xlabel = "Re u", ylabel = "Im u", aspect = DataAspect(),
              title = "The u-plane: wall of marginal stability {u : a_D/a ∈ ℝ}")
    lines!(ax, real.(closed), imag.(closed); color = :black, linewidth = 2)
    scatter!(ax, [2.0, -2.0], [0.0, 0.0]; color = :crimson, markersize = 14)
    text!(ax, 2.1, 0.12; text = "monopole  u = 2Λ²")
    text!(ax, -2.9, 0.12; text = "dyon  u = −2Λ²")
    text!(ax, 0.0, 0.0; text = "strong coupling\nmonopole + dyon", align = (:center, :center))
    text!(ax, 0.0, 2.3; text = "weak coupling\nW-boson + two dyon towers",
          align = (:center, :center))
    xlims!(ax, -4.2, 4.2); ylims!(ax, -3.2, 3.2)
    save_fig("sw_uplane.png", fig)
end

# Periods along the real axis and the phase whose lock to ℝ marks the wall
let us = range(2.2, 60.0; length = 300), ps = [sw_periods(sw, u) for u in us]
    fig = Figure(size = (900, 340))
    ax1 = Axis(fig[1, 1]; xlabel = "u", ylabel = "period",
               title = "Seiberg-Witten periods on the real axis")
    lines!(ax1, us, [abs(p.a) for p in ps]; linewidth = 2, label = "|a|  (electric)")
    lines!(ax1, us, [abs(p.a_D) for p in ps]; linewidth = 2, label = "|a_D|  (magnetic)")
    axislegend(ax1; position = :lt, framevisible = false)
    ax2 = Axis(fig[1, 2]; xlabel = "u", ylabel = "arg(a_D / a)",
               title = "period ratio: the phase that locks on the wall")
    lines!(ax2, us, [angle(p.a_D / p.a) for p in ps]; linewidth = 2)
    save_fig("sw_periods.png", fig)
end

# Central charges of the BPS states in both chambers. Rays are scaled by the heaviest
# state in each panel, so both the phase ordering and the mass ratios are visible.
let
    fig = Figure(size = (940, 460))
    # rays are drawn at unit length: what matters for stability is the phase of Z, and
    # the masses span two orders of magnitude in the weak-coupling tower
    circle_pts = [cis(t) for t in range(0, 2π; length = 200)]
    function draw_chamber!(ax, states; label_charges = false)
        lines!(ax, real.(circle_pts), imag.(circle_pts); color = (:gray, 0.35))
        for st in states
            Z = ComplexF64(central_charge(st))
            e = Z / abs(Z)
            vector = ExactWKB.omega(st) == -2
            lines!(ax, [0, real(e)], [0, imag(e)];
                   linewidth = vector ? 3.5 : 1.6,
                   color = vector ? :crimson : :steelblue)
            scatter!(ax, [real(e)], [imag(e)]; markersize = vector ? 12 : 7,
                     color = vector ? :crimson : :steelblue)
            label_charges && text!(ax, 1.12 * real(e), 1.12 * imag(e);
                                   text = string(charge(st)), fontsize = 14,
                                   align = (:center, :center))
        end
        xlims!(ax, -1.35, 1.35); ylims!(ax, -1.35, 1.35)
    end
    ax1 = Axis(fig[1, 1]; xlabel = "Re Z / |Z|", ylabel = "Im Z / |Z|",
               aspect = DataAspect(),
               title = "Strong coupling (u = 0.4i): monopole and dyon")
    draw_chamber!(ax1, su2_bps_states(sw, 0.4im); label_charges = true)
    ax2 = Axis(fig[1, 2]; xlabel = "Re Z / |Z|", ylabel = "Im Z / |Z|",
               aspect = DataAspect(),
               title = "Weak coupling (u = 25): towers accumulate on the W-boson")
    draw_chamber!(ax2, su2_bps_states(sw, 25.0 + 0im; tower = 8))
    text!(ax2, -0.02, 0.13; text = "W-boson (Ω = −2)", color = :crimson, fontsize = 13,
          align = (:left, :center))
    save_fig("sw_bps_charges.png", fig)
end

# The Nekrasov-Shatashvili quantum correction to the electric period
let us = range(6.0, 60.0; length = 120)
    a2 = [real(Resurgence.coefficients(quantum_sw_periods(sw, u; order = 1).a)[3]) for u in us]
    fig = Figure(size = (660, 420))
    ax = Axis(fig[1, 1]; xlabel = "u", ylabel = "−a₂(u)", yscale = log10,
              title = "Quantum correction to the electric period and its large-u law")
    lines!(ax, us, -a2; linewidth = 2, label = "−a₂ (Picard-Fuchs reduced operator)")
    lines!(ax, us, [1 / (4u^2.5) for u in us]; color = :orangered, linestyle = :dash,
           linewidth = 2, label = "Λ⁴ / (4 u^{5/2})")
    axislegend(ax; position = :rt, framevisible = false)
    save_fig("sw_quantum.png", fig)
end

println("\nAll figures written to ", ASSETS)
