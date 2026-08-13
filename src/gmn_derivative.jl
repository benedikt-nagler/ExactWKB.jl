# Moduli derivatives of the finite-radius GMN solution, by implicit differentiation of the
# fixed point instead of finite differences of the solver. `metric_point` used nine
# independent `solve_gmn` calls at h = 1e-3 to build ∂_μ log 𝒳; with a basis period
# derivative dZ/du in hand the eight neighbours collapse into one linear solve, and the
# metric's accuracy stops being set by h.
#
# The discrete fixed point of `solve_gmn` carries the moduli in exactly two places, since
# `metric_point` shares one rapidity grid across the evaluation:
#
#   1. the semiflat source −2πR|Z_a| cosh θ_k + iθ_{γ_a},
#   2. the rays e_a = −Z_a/|Z_a| inside the Toeplitz kernel.
#
# Both are elementary. |Z| gives ∂_μ|Z_a| = Re(Z̄_a ∂_μZ_a)/|Z_a|; the ray has unit modulus,
# so ∂_μ e_a = iλ_a e_a with λ_a = Im(∂_μZ_a/Z_a) real, and with p = ρe_b, q = e_a,
#
#   ∂_μ A = 2i(λ_a − λ_b)·pq/(p−q)²,     ∂_μ B = −2i(λ_a − λ_b)·pq/(p+q)²
#
# which vanishes when every ray turns together, as a rigid rotation of the ζ-plane must.
# The linearization is then (I − M)δg = b_μ with M one existing Picard sweep carrying
# ∂_g log1p(−σe^g) in place of log1p(−σe^g). It is solved matrix-free by Picard on the
# *linear* equation, all 2 + rank right-hand sides in one loop: the operator M is the
# derivative of a contraction that converged, so the linear iteration contracts at the same
# rate. Second derivatives are the same operator with a different right-hand side, which is
# the argument for the implicit function theorem over an AD Hessian here.
#
# ── Convention ledger (each item pinned by an oracle in test/test_gmn_derivative.jl) ────
#
# 1. Moduli order (Re u, Im u, θ₁, …, θ_r), the ordering `hyperkahler.jl`'s metric block
#    already uses. Every derivative here is indexed μ = 1 … 2 + r in that order.
# 2. `dZdu` is the **holomorphic** u-derivative of the *basis* central charges, so the two
#    real directions are ∂_{Re u}Z = Z′ and ∂_{Im u}Z = iZ′. State central charges follow
#    by linearity, Z_a = Σ_i γ_a[i] Z_{e_i}, which `gmn_torus` has already checked. Pinned
#    by the Ooguri-Vafa oracle, where the two real directions are compared separately
#    against a quadrature independent of `solve_gmn`.
# 3. The BPS data (charges, Ω, σ, ⟨,⟩) is held **constant**: it is locally constant integer
#    data away from a wall of marginal stability. On the wall the derivative is undefined,
#    and so is the finite difference it replaces.

# ── the linearized nonlinearity ───────────────────────────────────────────────────────

# ∂_g of `_gmn_lplus` and `_gmn_lminus`
_dgmn_lplus(g, σ) = (x = σ * exp(g); -x / (1 - x))
_dgmn_lminus(g, σ) = (x = σ * exp(-g); x / (1 - x))

# ∂A and ∂B stripped of their scalar (λ_a − λ_b): the Toeplitz partners of `_gmn_kernels`
function _gmn_ray_kernels(t::GMNTorus{F}, np::Integer, dθ::F) where {F}
    n = n_states(t)
    ρs = [exp(d * dθ) for d in -(np - 1):(np - 1)]
    GA = Dict{Tuple{Int,Int},Vector{Complex{F}}}()
    GB = Dict{Tuple{Int,Int},Vector{Complex{F}}}()
    for a in 1:n, b in 1:n
        t.state_pairing[a, b] == 0 && continue
        q = t.rays[a]
        GA[(a, b)] = Complex{F}[(p = ρ * t.rays[b]; 2 * im * p * q / (p - q)^2) for ρ in ρs]
        GB[(a, b)] = Complex{F}[(p = ρ * t.rays[b]; -2 * im * p * q / (p + q)^2) for ρ in ρs]
    end
    GA, GB
end

"""
    GMNDerivative{F}

The moduli derivatives of a [`GMNSolution`](@ref): `dgplus[k,a,μ]` and `dgminus[k,a,μ]`
hold ``∂_μ \\log 𝒳_{γ_a}`` at the `k`-th grid point of the ray of state `a`, for the
moduli ``μ = (\\mathrm{Re}\\,u, \\mathrm{Im}\\,u, θ_1, …, θ_r)``. Build with
[`solve_gmn_derivative`](@ref) and evaluate off the grid with
[`log_xi_derivative`](@ref). `n_iterations` and `residual` report the linear solve.
"""
struct GMNDerivative{F}
    solution::GMNSolution{F}
    dZdu::Vector{Complex{F}}
    lambda::Matrix{F}
    dgplus::Array{Complex{F},3}
    dgminus::Array{Complex{F},3}
    iterations::Int
    residual::F
end

n_iterations(gd::GMNDerivative) = gd.iterations
residual(gd::GMNDerivative) = gd.residual
n_states(gd::GMNDerivative) = n_states(gd.solution)
n_charges(gd::GMNDerivative) = n_charges(gd.solution.torus)

"""
    n_parameters(gd::GMNDerivative) -> Int

The number of real moduli differentiated, ``2 + r`` for a rank-``r`` charge lattice.
"""
n_parameters(gd::GMNDerivative) = size(gd.dgplus, 3)

"""
    solve_gmn_derivative(sol::GMNSolution, dZdu::AbstractVector; tol = nothing,
                         maxiter = 500, relax = 1) -> GMNDerivative

Differentiate a converged GMN solution with respect to the moduli
``(\\mathrm{Re}\\,u, \\mathrm{Im}\\,u, θ_1, …, θ_r)``, by linearizing the integral
equation rather than by re-solving it at shifted moduli. `dZdu` gives ``dZ/du`` for the
`r` charges of the lattice basis, in the basis `sol`'s charges are written in;
[`sw_period_derivatives`](@ref) and [`ooguri_vafa_period_derivatives`](@ref) supply it in
closed form.

The linearized equation ``(I - M)δg = b_μ`` is solved by iteration on the same Toeplitz
kernel the solution used, with all ``2 + r`` right-hand sides in one loop. `M` is the
derivative of a map that converged, so the iteration contracts at the rate `sol` did;
`tol` is a **relative** sup-norm change. Throws [`TBAError`](@ref) if it has not
converged within `maxiter` sweeps.

The BPS spectrum is held fixed, so the derivative is the one-sided limit on a wall of
marginal stability, where the spectrum jumps. Convention: `dZdu` is holomorphic in `u`,
see the ledger at the top of `src/gmn_derivative.jl`.
"""
function solve_gmn_derivative(sol::GMNSolution{F}, dZdu::AbstractVector;
                              tol = nothing, maxiter::Integer = 500,
                              relax::Real = 1) where {F}
    0 < relax ≤ 1 || throw(Resurgence.InvalidArgument(
        "relax must be in (0, 1], got $relax"))
    t = sol.torus
    r = n_charges(t)
    length(dZdu) == r || throw(Resurgence.InvalidArgument(
        "expected $r basis period derivatives dZ/du, got $(length(dZdu))"))
    n = n_states(t)
    npar = 2 + r
    tol = tol === nothing ? sqrt(eps(F)) * F(1 // 100) : F(tol)
    grid, gp, gm = sol.grid, sol.gplus, sol.gminus
    np = length(grid)
    w = _trapezoid(grid)
    m = abs.(t.Z)
    dZb = Complex{F}[Complex{F}(z) for z in dZdu]
    dZ = Complex{F}[sum(t.charges[a][i] * dZb[i] for i in 1:r) for a in 1:n]

    # the moduli enter through the ray angle (λ), the mass (dm) and the torus angle (dang)
    λ = zeros(F, n, npar)
    dm = zeros(F, n, npar)
    dang = zeros(F, n, npar)
    for a in 1:n
        for (μ, dir) in ((1, one(Complex{F})), (2, Complex{F}(im)))
            d = dir * dZ[a]
            λ[a, μ] = imag(d / t.Z[a])
            dm[a, μ] = real(conj(t.Z[a]) * d) / m[a]
        end
        for i in 1:r
            dang[a, 2 + i] = t.charges[a][i]
        end
    end

    dθ = grid[2] - grid[1]
    KA, KB = _gmn_kernels(t, np, dθ)
    GA, GB = _gmn_ray_kernels(t, np, dθ)

    # the source: ∂_μ(semiflat) plus the motion of the rays through the kernel
    Lp = Matrix{Complex{F}}(undef, np, n)
    Lm = Matrix{Complex{F}}(undef, np, n)
    for b in 1:n
        σ = t.sigma[b]
        @views @. Lp[:, b] = w * _gmn_lplus(gp[:, b], σ)
        @views @. Lm[:, b] = w * _gmn_lminus(gm[:, b], σ)
    end
    bp = Array{Complex{F},3}(undef, np, n, npar)
    bm = Array{Complex{F},3}(undef, np, n, npar)
    for μ in 1:npar, a in 1:n, k in 1:np
        s = 2 * F(π) * t.R * cosh(grid[k]) * dm[a, μ]
        bp[k, a, μ] = -s + im * dang[a, μ]
        bm[k, a, μ] = +s + im * dang[a, μ]
    end
    for a in 1:n, b in 1:n
        c = t.state_pairing[a, b]
        c == 0 && continue
        coef = t.omega[b] * c / (4 * F(π) * im)
        A, B = GA[(a, b)], GB[(a, b)]
        for μ in 1:npar
            dλ = λ[a, μ] - λ[b, μ]
            iszero(dλ) && continue
            cf = coef * dλ
            for k in 1:np
                accp = zero(Complex{F})
                accm = zero(Complex{F})
                @inbounds for l in 1:np
                    d = k - l + np
                    accp += A[d] * Lp[l, b] - B[d] * Lm[l, b]
                    accm += B[d] * Lp[l, b] - A[d] * Lm[l, b]
                end
                bp[k, a, μ] += cf * accp
                bm[k, a, μ] += cf * accm
            end
        end
    end

    # (I − M)δ = b by iteration; M is one sweep with the differentiated nonlinearity
    δp = copy(bp)
    δm = copy(bm)
    newp = similar(δp)
    newm = similar(δm)
    Pp = similar(δp)
    Pm = similar(δm)
    dLp = Matrix{Complex{F}}(undef, np, n)
    dLm = Matrix{Complex{F}}(undef, np, n)
    for b in 1:n
        σ = t.sigma[b]
        @views @. dLp[:, b] = w * _dgmn_lplus(gp[:, b], σ)
        @views @. dLm[:, b] = w * _dgmn_lminus(gm[:, b], σ)
    end
    res = F(Inf)
    iters = 0
    for it in 1:maxiter
        iters = it
        for μ in 1:npar, b in 1:n, l in 1:np
            Pp[l, b, μ] = dLp[l, b] * δp[l, b, μ]
            Pm[l, b, μ] = dLm[l, b] * δm[l, b, μ]
        end
        copyto!(newp, bp)
        copyto!(newm, bm)
        for a in 1:n, b in 1:n
            c = t.state_pairing[a, b]
            c == 0 && continue
            coef = t.omega[b] * c / (4 * F(π) * im)
            A, B = KA[(a, b)], KB[(a, b)]
            for μ in 1:npar, k in 1:np
                accp = zero(Complex{F})
                accm = zero(Complex{F})
                @inbounds for l in 1:np
                    d = k - l + np
                    accp += A[d] * Pp[l, b, μ] - B[d] * Pm[l, b, μ]
                    accm += B[d] * Pp[l, b, μ] - A[d] * Pm[l, b, μ]
                end
                newp[k, a, μ] += coef * accp
                newm[k, a, μ] += coef * accm
            end
        end
        scale = max(maximum(abs, δp), maximum(abs, δm), one(F))
        res = max(maximum(abs, newp - δp), maximum(abs, newm - δm)) / scale
        @. δp = δp + relax * (newp - δp)
        @. δm = δm + relax * (newm - δm)
        res < tol && return GMNDerivative{F}(sol, dZb, λ, δp, δm, iters, res)
    end
    throw(TBAError("the linearized GMN iteration did not converge: relative sup-norm " *
                   "change $res after $maxiter sweeps (tol $tol) - pass a larger " *
                   "maxiter, or relax < 1"))
end

"""
    log_xi_derivative(gd::GMNDerivative, γ::AbstractVector{<:Integer}, ζ) -> Vector

The moduli gradient ``∂_μ \\log 𝒳_γ(ζ)`` of the Darboux coordinate at **fixed** ``ζ``,
one entry per modulus in the order
``(\\mathrm{Re}\\,u, \\mathrm{Im}\\,u, θ_1, …, θ_r)``. This is the exact form of the
central difference [`metric_point`](@ref) takes without a `dZdu`.

The logarithm is differentiated, not ``𝒳_γ`` itself, which is what stays finite at large
``|ζ|`` and what the metric consumes.

Throws [`TBAError`](@ref) if `ζ` lies on a BPS ray. There ``𝒳_γ`` jumps by the
Kontsevich-Soibelman factor and the one-sided derivatives differ.
"""
function log_xi_derivative(gd::GMNDerivative{F}, γ::AbstractVector{<:Integer},
                           ζ::Number) where {F}
    sol = gd.solution
    t = sol.torus
    r = n_charges(t)
    length(γ) == r || throw(Resurgence.InvalidArgument(
        "charge $γ does not have the lattice rank $r"))
    iszero(ζ) && throw(TBAError("𝒳_γ has an essential singularity at ζ = 0"))
    npar = 2 + r
    n = n_states(t)
    grid, gp, gm = sol.grid, sol.gplus, sol.gminus
    ts = exp.(-grid)
    w = _trapezoid(grid)

    out = zeros(Complex{F}, npar)
    dZγ = sum(γ[i] * gd.dZdu[i] for i in 1:r)
    for (μ, dir) in ((1, one(Complex{F})), (2, Complex{F}(im)))
        d = dir * dZγ
        out[μ] = F(π) * t.R * (d / ζ + conj(d) * ζ)
    end
    for i in 1:r
        out[2 + i] = im * γ[i]
    end
    for b in 1:n
        c = transpose(collect(γ)) * t.pairing * t.charges[b]
        c == 0 && continue
        σb = t.sigma[b]
        coef = t.omega[b] * c / (4 * F(π) * im)
        for ε in (1, -1)
            ej = ε * t.rays[b]
            abs(angle(ζ / ej)) < sqrt(eps(F)) && throw(TBAError(
                "ζ = $ζ lies on the BPS ray of state $b, where 𝒳_γ jumps and the " *
                "moduli derivative is one-sided - move off the ray"))
            col = ε == 1 ? view(gp, :, b) : view(gm, :, b)
            dcol = ε == 1 ? view(gd.dgplus, :, b, :) : view(gd.dgminus, :, b, :)
            @inbounds for k in eachindex(ts)
                ζp = ts[k] * ej
                K = _tba_kernel(ζ, ζp)
                dK = -2 * ζ * ζp / (ζp - ζ)^2
                L = ε == 1 ? _gmn_lplus(col[k], σb) : _gmn_lminus(col[k], σb)
                dL = ε == 1 ? _dgmn_lplus(col[k], σb) : _dgmn_lminus(col[k], σb)
                for μ in 1:npar
                    out[μ] += ε * coef * w[k] *
                              (dK * im * gd.lambda[b, μ] * L + K * dL * dcol[k, μ])
                end
            end
        end
    end
    out
end
