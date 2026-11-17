# TBA flagship demonstration: exact quantum periods from BPS data alone. The
# conformal-limit Gaiotto–Moore–Neitzke / Ito–Mariño–Shu integral equations are solved
# from nothing but the spectrum (Z_γ, Ω(γ), ⟨γ,γ'⟩) and compared against the entirely
# independent series pipeline (WKB recursion → Borel–Padé median sums) on the cubic
# and quartic potentials - plus the constant-Y-system plateau (the golden ratio) and
# the lateral jump reproducing the DDP wall-crossing formula.
#
# Run:  julia --project=test examples/tba_demo.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "test"))

using ExactWKB
import Resurgence

function borel_side(prob, w, Z, ħ)
    tps = simple_turning_points(prob)
    for i in 1:(length(tps) - 1), j in (i + 1):length(tps)
        c = try
            encircling_contour(tps[i], tps[j]; margin = 0.4)
        catch
            continue
        end
        cp = try
            period_integral(prob, c)
        catch
            continue
        end
        s = abs(cp - Z) < 1e-6 * abs(Z) ? 1 : abs(cp + Z) < 1e-6 * abs(Z) ? -1 : 0
        s == 0 && continue
        return s * ExactWKB._log_voros(voros_symbol(w, c), ħ)
    end
    nothing
end

for (name, coeffs) in (("cubic  Q = z³ − z", [0.0, -1.0, 0.0, 1.0]),
                       ("quartic Q = z⁴ − 1", [-1.0, 0.0, 0.0, 0.0, 1.0]))
    println("── ", name, " ", "─"^(50 - length(name)))
    prob = SchrodingerProblem(coeffs)
    sp = bps_spectrum(prob)
    println(sp)
    sol = solve_tba(sp)
    println(sol)
    w = wkb_expansion(prob; order = 12)
    println("\n  log X_γ(ħ) at ħ = |Z_γ|/6 - TBA (no series) vs Borel–Padé median:")
    for a in 1:n_states(sp)
        Z = sol.system.Z[a]
        ħ = abs(Z) / 6
        lt = ExactWKB._log_x(sol, a, ħ, 0)
        lb = borel_side(prob, w, Z, ħ)
        lb === nothing && continue
        rel = abs(lt - lb) / abs(lb)
        println("    γ_", a, ":  TBA = ", round(lt; sigdigits = 10),
                "\n         Borel = ", round(lb; sigdigits = 10),
                "   (rel. diff ", round(rel; sigdigits = 2), ")")
        @assert rel < 1e-3
    end
    println()
end

# the IR plateau of the cubic solves Zamolodchikov's constant A₂ Y-system: X² = 1 + X
cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
sol = solve_tba(bps_spectrum(cubic))
k = findfirst(θ -> θ > sol.grid[1] + 12, sol.grid)
x = real(exp(sol.gplus[k, 1]))
println("IR plateau of the cubic: X(θ → −∞) = ", round(x; sigdigits = 8),
        "  vs golden ratio ", round((1 + sqrt(5)) / 2; sigdigits = 8))
@assert abs(x - (1 + sqrt(5)) / 2) < 5e-3

# the lateral jump across the θ_c = 0 ray is the DDP / Kontsevich–Soibelman factor
sys = sol.system
i = argmax([abs(imag(z)) for z in sys.Z])
j = 3 - i
ħ = abs(sys.Z[j]) / 6
jump = ExactWKB._log_x(sol, i, ħ, 1) - ExactWKB._log_x(sol, i, ħ, -1)
κ = jump / log1p(voros_value(sol, j, ħ))
println("Stokes constant from the TBA lateral jump: κ = ", round(κ; sigdigits = 6),
        "  (exact integer ⟨γ,γ'⟩Ω = ", sys.pairing[i, j] * sys.omega[j], ")")
@assert abs(κ - sys.pairing[i, j] * sys.omega[j]) < 1e-2
println("\nAll TBA demonstration assertions passed.")
