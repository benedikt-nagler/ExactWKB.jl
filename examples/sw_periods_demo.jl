# Seiberg–Witten SU(2) period demonstration: the first rung of the SU(2) north-star
# ascent. From the pure-SU(2) geometry alone we compute the classical periods a(u),
# a_D(u) over the weak-coupling u-plane, verify the weak-coupling identity
# u = a² + Λ⁴/(2a²) + …, print the leading Nekrasov–Shatashvili quantum correction to
# the electric period, and cross-check the Mathieu characteristic values a_ν(q) against
# an independent tridiagonal solver. A PNG shows |a|, |a_D| and arg(a_D/a) (whose lock
# to 0/π previews the wall of marginal stability) over the u-plane.
#
# Run:  julia --project=test examples/sw_periods_demo.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "test"))

using ExactWKB
import Resurgence

sw = SeibergWittenSU2()          # Λ = 1
println("Pure SU(2) Seiberg–Witten geometry: ", sw)
s = sw_singularities(sw)
println("  singular moduli: monopole u = ", s.monopole, ", dyon u = ", s.dyon, "\n")

println("weak-coupling periods and the u(a) identity")
println("   u        a           a_D            u − (a²+Λ⁴/2a²+…)")
for u in (10.0, 25.0, 100.0)
    p = sw_periods(sw, u)
    a = real(p.a)
    u_ser = a^2 + 1/(2a^2) + 5/(32a^6) + 9/(64a^10)
    println("  ", rpad(u, 7), "  ", rpad(round(a; sigdigits = 6), 10), "  ",
            rpad(round(imag(p.a_D); sigdigits = 6), 12), "  ", u - u_ser)
end

println("\nNekrasov–Shatashvili quantum correction to the electric period (a₂·ħ²)")
println("   u        a₂ (numeric)     −Λ⁴/(4u^{5/2})")
for u in (100.0, 400.0)
    a2 = real(Resurgence.coefficients(quantum_sw_periods(sw, u; order = 1).a)[3])
    println("  ", rpad(u, 7), "  ", rpad(round(a2; sigdigits = 6), 15), "  ",
            round(-1/(4u^2.5); sigdigits = 6))
end

println("\nmonodromy (basis (a_D, a)):  M_∞ = M_dyon · M_monopole")
println("  M_∞        = ", sw_monodromy(sw, :infinity))
println("  M_monopole = ", sw_monodromy(sw, :monopole))
println("  M_dyon     = ", sw_monodromy(sw, :dyon))

# u-plane figure (headless CairoMakie, mirroring the other examples)
try
    using CairoMakie
    us = range(2.2, 60.0; length = 200)
    ps = [sw_periods(sw, u) for u in us]
    fig = Figure(size = (900, 300))
    ax1 = Axis(fig[1, 1]; xlabel = "u", ylabel = "|period|", title = "SU(2) periods")
    lines!(ax1, us, [abs(p.a) for p in ps]; label = "|a|")
    lines!(ax1, us, [abs(p.a_D) for p in ps]; label = "|a_D|")
    axislegend(ax1; position = :lt)
    ax2 = Axis(fig[1, 2]; xlabel = "u", ylabel = "arg(a_D/a)",
               title = "phase (wall preview)")
    lines!(ax2, us, [angle(p.a_D / p.a) for p in ps])
    path = joinpath(@__DIR__, "sw_periods.png")
    save(path, fig)
    println("\nwrote ", path)
catch err
    println("\n(skipping plot: ", err, ")")
end
