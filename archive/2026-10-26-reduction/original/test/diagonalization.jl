# Test-side eigenvalue oracle (NOT package code): dense diagonalization of the
# truncated Hamiltonian H = −ħ²d² + V(z) in the harmonic-oscillator basis of
# H₀ = −ħ²d² + ω²z². LinearAlgebra stdlib only. Float64 LAPACK eigenvalues seed a
# BigFloat Rayleigh-quotient (inverse) iteration - generic `eigen` has no BigFloat
# dense path in the stdlib, but generic `lu` does, and RQI converges cubically from
# a LAPACK-accurate seed. Target: 20+ digits for the quartic / double well.
#
# Ladder conventions (matching ħ²ψ″ = (V − E)ψ, i.e. −ħ²ψ″ + Vψ = Eψ):
#   z = √(ħ/(2ω))·(a + a†),   −ħ²d² = p²,  p = i√(ħω/2)·(a† − a),
# so H₀ has eigenvalues 2ħω(n + ½) (the ω = 1 harmonic gives E_n = 2ħ(n + ½)).

using LinearAlgebra

# H as a dense symmetric matrix in the first N basis states, element type F
function _hamiltonian_matrix(v_coeffs, ħ, ω, N, ::Type{F}) where {F}
    a = zeros(F, N, N)
    for k in 1:(N - 1)
        a[k, k + 1] = sqrt(F(k))          # a|n⟩ = √n |n−1⟩ (1-based rows)
    end
    zmat = sqrt(F(ħ) / (2 * F(ω))) .* (a .+ a')
    p2 = -(F(ħ) * F(ω) / 2) .* ((a' .- a) * (a' .- a))
    # V(zmat) by Horner on the ascending coefficients
    V = F(last(v_coeffs)) .* Matrix{F}(I, N, N)
    for k in (length(v_coeffs) - 1):-1:1
        V = zmat * V + F(v_coeffs[k]) .* Matrix{F}(I, N, N)
    end
    H = p2 .+ V
    (H .+ H') ./ 2
end

# Rayleigh-quotient iteration in F, seeded by (λ0, x0)
function _rqi(H::Matrix{F}, λ0, x0; iters = 8) where {F}
    λ = F(λ0)
    x = F.(x0)
    x ./= norm(x)
    for _ in 1:iters
        y = try
            lu(H - λ * I) \ x
        catch e
            e isa SingularException ? x : rethrow()   # λ exact: converged
        end
        x = y ./ norm(y)
        λ = dot(x, H * x)
    end
    λ
end

"""
    diagonalization_eigenvalues(prob, ħ; nev = 4, N = 120, omega = 1.0,
                                refine = true) -> Vector

The lowest `nev` eigenvalues of `−ħ²ψ″ + Vψ = Eψ` by dense diagonalization of the
`N`-state harmonic-basis truncation. Float64 LAPACK first; with `refine = true`
each level is polished by BigFloat Rayleigh-quotient iteration at the caller's
`setprecision`. Choose `N` large enough that the returned values are stable - the
truncation error, not arithmetic, is the accuracy limit.
"""
function diagonalization_eigenvalues(prob, ħ; nev::Integer = 4, N::Integer = 120,
                                     omega::Real = 1.0, refine::Bool = true)
    v = Float64.(prob.v_coeffs)
    H64 = _hamiltonian_matrix(v, Float64(ħ), Float64(omega), N, Float64)
    ev = eigen(Symmetric(H64))
    λs = ev.values[1:nev]
    refine || return λs
    vb = BigFloat.(prob.v_coeffs)
    Hb = _hamiltonian_matrix(vb, BigFloat(ħ), BigFloat(omega), N, BigFloat)
    [_rqi(Hb, λs[k], ev.vectors[:, k]) for k in 1:nev]
end
