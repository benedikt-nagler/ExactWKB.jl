# SU(2) BPS spectrum + wall-crossing demonstration (rung 2 of the SU(2) ascent). The two BPS
# chambers of pure SU(2) - strong coupling (monopole + dyon) and weak coupling (W-boson +
# infinite dyon tower) - are printed with their charges, DT invariants Ω and masses (from the
# rung-1 Seiberg–Witten periods), the affine-Kronecker BPS quiver is exercised against
# ClusterAlgebras, the Kontsevich–Soibelman wall-crossing identity is verified at increasing
# degree, and the dyon tower's phase accumulation onto the W-boson is tabulated.
#
# Run:  julia --project=test examples/su2_bps_demo.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "test"))

using ExactWKB
const CA = ExactWKB.ClusterAlgebras

sw = SeibergWittenSU2()
u = 25.0
println("Pure SU(2) BPS spectrum at u = ", u, "  (", sw, ")\n")

function show_states(label, states)
    println(label)
    println("   charge (n_m,n_e)     Ω     |Z|          arg Z")
    for s in states
        Z = ExactWKB.central_charge(s)
        println("   ", rpad(Tuple(ExactWKB.charge(s)), 18), " ", rpad(ExactWKB.omega(s), 5),
                " ", rpad(round(ExactWKB.mass(s); sigdigits = 6), 12),
                " ", round(angle(Z); sigdigits = 4))
    end
    println()
end

show_states("strong coupling (inside the wall):", su2_bps_states(sw, u; chamber = :strong))
show_states("weak coupling (outside the wall):", su2_bps_states(sw, u; chamber = :weak, tower = 3))

q = su2_bps_quiver()
println("BPS quiver = Kronecker: is_affine=", CA.is_affine_type(q),
        "  is_finite=", CA.is_finite_type(q),
        "  is_mutation_finite=", CA.is_mutation_finite(q), "\n")

println("Kontsevich–Soibelman wall-crossing (strong ⇔ weak), by truncation degree:")
for d in 2:6
    r = verify_su2_wall_crossing(; degree = d)
    println("   degree ", d, ":  closed = ", r.closed, "   (max residual ", r.max_residual, ")")
end

println("\ndyon-tower phase accumulation onto arg Z_W = ", round(angle(central_charge(sw, u, (0, 2))); sigdigits = 4))
println("   n     arg Z(1,2n)      arg Z(-1,2n+2)")
for n in (1, 2, 4, 8, 16, 32)
    pp = angle(central_charge(sw, u, (1, 2n)))
    pi_ = angle(central_charge(sw, u, (-1, 2n + 2)))
    println("   ", rpad(n, 5), " ", rpad(round(pp; sigdigits = 4), 15), " ", round(pi_; sigdigits = 4))
end
