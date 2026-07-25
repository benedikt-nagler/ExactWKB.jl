# Seiberg–Witten SU(2): the wall of marginal stability on the u-plane. From the
# closed-form periods (complex u, Carlson elliptic integrals) we trace the wall
# {u : a_D/a ∈ ℝ} through the singular points ±2Λ², verify the two-sided BPS chamber
# structure it separates (monopole + dyon inside; W-boson + dyon towers outside),
# and check the Picard–Fuchs monodromy loops against the exact integer matrices.
# A PNG shows the wall, the singularities, and chamber-labeled sample points.
#
# Run:  julia --project=test examples/su2_wall_demo.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "test"))

using ExactWKB

sw = SeibergWittenSU2()          # Λ = 1
s = sw_singularities(sw)
println("Pure SU(2): singular moduli u = ", s.monopole, " (monopole), ", s.dyon, " (dyon)\n")

wall = ms_wall(sw; n = 64)
crossing = wall[2 + 64 ÷ 2]      # nearest sampled point to the imaginary axis
println("wall of marginal stability {u : a_D/a ∈ ℝ}: ", length(wall), " points,")
println("  through ", wall[1], " and ", wall[66], ", crossing near iℝ at ", crossing)
r = sw_periods(sw, crossing).a_D / sw_periods(sw, crossing).a
println("  there a_D/a = ", r, "  (real, in (−2, 0))\n")

println("chambers and BPS spectra (chamber = :auto)")
for u in (0.4im, 25.0 + 0im)
    ch = sw_chamber(sw, u)
    states = su2_bps_states(sw, u; tower = 3)
    println("  u = ", u, "  →  ", ch, " chamber, ", length(states), " states:")
    for st in states
        println("      charge ", ExactWKB.charge(st), "   Ω = ", ExactWKB.omega(st),
                "   |Z| = ", round(ExactWKB.mass(st); sigdigits = 4))
    end
end

println("\nPicard–Fuchs monodromy loops (transport vs exact integer matrices)")
base_u = 3.0 + 0im
circle(c, ρ, K) = [c + ρ * cis(2π * k / K) for k in 0:K]
function loop_matrix(path)
    b = continue_periods(sw, [base_u]); S = [b.a_D b.a; b.da_D b.da]
    f = continue_periods(sw, path);     S2 = [f.a_D f.a; f.da_D f.da]
    round.(Int, real.(transpose(S \ S2)))
end
am = [base_u, 2.6 + 0im]
ad = [base_u, 2.6im, -1.4 + 0im]
loops = (infinity = circle(0.0 + 0im, 3.0, 48),
         monopole = vcat(am, circle(2.0 + 0im, 0.6, 48)[2:end], reverse(am)),
         dyon = vcat(ad, circle(-2.0 + 0im, 0.6, 48)[2:end], reverse(ad)))
for pt in (:infinity, :monopole, :dyon)
    M = loop_matrix(loops[pt])
    println("  ", rpad(pt, 9), " transport = ", M, "   exact = ", sw_monodromy(sw, pt),
            "   match: ", M == sw_monodromy(sw, pt))
end

# u-plane figure (headless CairoMakie, mirroring the other examples)
try
    using CairoMakie
    fig = Figure(size = (600, 500))
    ax = Axis(fig[1, 1]; xlabel = "Re u", ylabel = "Im u",
              title = "SU(2) wall of marginal stability  {u : a_D/a ∈ ℝ}")
    closed = vcat(wall, wall[1:1])
    lines!(ax, real.(closed), imag.(closed); color = :black)
    scatter!(ax, [2.0, -2.0], [0.0, 0.0]; color = :red, markersize = 12)
    text!(ax, 2.05, 0.1; text = "monopole")
    text!(ax, -1.95, 0.1; text = "dyon")
    text!(ax, -0.45, 0.0; text = "strong:\nmonopole + dyon", align = (:center, :center))
    text!(ax, 0.0, 1.55; text = "weak: W-boson + dyon towers", align = (:center, :center))
    path = joinpath(@__DIR__, "su2_wall.png")
    save(path, fig)
    println("\nwrote ", path)
catch err
    println("\n(skipping plot: ", err, ")")
end
