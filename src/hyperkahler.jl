# Hyperkähler metrics by the Gaiotto–Moore–Neitzke twistor construction [GMN10]. This is
# the finite-radius companion of `tba.jl`: the same integral equation, but with the radius
# R of the circle the 4d theory is compactified on, and with the torus angles θ_γ restored.
# The conformal limit implemented in `tba.jl` drops both and gives quantum periods; keeping
# both gives Darboux coordinates 𝒳_γ(ζ) on the twistor space of the moduli space M of the
# theory on ℝ³ × S¹, and from them the hyperkähler metric on M - semiflat plus the full
# instanton series, resummed.
#
# Everything here consumes **lattice-level** BPS data only (charges, Ω, Z, ⟨,⟩): a
# `Vector{BPSState}` plus the Gram matrix of the basis those charges are written in. That
# is what makes it work for both producers in this package - `bps_spectrum` (charges in a
# `ChargeBasis`, Gram `signed_pairing`) and `su2_bps_states` (physical charges (n_m, n_e),
# Gram [0 1; −1 0]).
#
# ── Convention ledger (each item pinned by an oracle in test/test_hyperkahler.jl) ──────
#
# 1. The equation [GMN10] (5.13). For every charge γ in the lattice:
#      log 𝒳_γ(ζ) = πR Z_γ/ζ + iθ_γ + πR Z̄_γ ζ
#          + (1/4πi) Σ_{γ'} Ω(γ') ⟨γ,γ'⟩ ∫_{ℓ_{γ'}} dζ'/ζ' (ζ'+ζ)/(ζ'−ζ)
#                                          · log(1 − σ(γ') 𝒳_{γ'}(ζ')),
#    sum over ±(spectrum charges), ray ℓ_{γ'} = {ζ : Z_{γ'}/ζ ∈ ℝ_{<0}}. The kernel and
#    the ±γ' folding are literally `tba.jl`'s (`_tba_kernel`, `_pair_sum`) - only the
#    inhomogeneous term differs. Pinned by the R → 0 oracle against `solve_tba`.
# 2. Semiflat normalization: πR Z/ζ + iθ_γ + πR Z̄ζ (GMN's, with 2πR absent from the
#    angle). On ℓ_γ with ζ = e^{−θ}e_γ, e_γ = −Z_γ/|Z_γ|, this is −2πR|Z_γ| cosh θ + iθ_γ,
#    so the rapidity window is **symmetric** [−θmax, θmax] (the conformal limit's window
#    is one-sided). Pinned by the Ω = 0 semiflat oracle and by the Ooguri–Vafa Bessel
#    oracle, which fixes the 2πR in the exponent K_0(2πR n |Z|).
# 3. Quadratic refinement σ(γ) ∈ {±1}, entering as log(1 − σ𝒳). Default σ ≡ −1, i.e.
#    log(1 + 𝒳) - exactly the conformal-limit convention of the `tba.jl` ledger, which is
#    what makes item 1's R → 0 oracle a like-for-like comparison. σ(−γ) = σ(γ).
# 4. Torus angles are linear in the charge: θ_γ = Σ_i γ_i θ_i with θ ∈ ℝ^rank given in
#    the same basis as the charges. Pinned by the semiflat oracle (𝒳_γ 𝒳_{γ'} = 𝒳_{γ+γ'}).
# 5. Charge pairing: whatever Gram matrix the caller passes, used as ⟨γ_a,γ_b⟩ = γ_aᵀ P γ_b.
#    For `su2_bps_states` that is the **physical Dirac pairing** [0 1; −1 0],
#    ⟨(m,e),(m',e')⟩ = m e' − e m' - which is `sw_bps.jl`'s ledger form and *equals* the
#    factor-2 form `_ks_skew` uses, because the latter is written in the abstract Kronecker
#    basis γ₁ = (1,0), γ₂ = (−1,2) where ⟨γ₁,γ₂⟩ = 2. Pinned by `su2_pairing`'s test.
# 6. Holomorphic symplectic form ϖ(ζ) = (1/4π²R) dlog𝒳_{γ₁} ∧ dlog𝒳_{γ₂} for a
#    **unimodular symplectic** basis ⟨γ₁,γ₂⟩ = 1, and the twistor theorem says it is a
#    Laurent polynomial ϖ = (−i/2ζ)ω₊ + ω₃ − (i/2)ζ ω₋. Both the prefactor and the
#    expansion are pinned by the semiflat limit, where ϖ is computed in closed form.
# 7. Metric from the hyperkähler triple: J_i = g⁻¹ω_i with J₁J₂ = ±J₃ forces
#    g = ω₁ ω₃⁻¹ ω₂ up to that sign, and **g = ω₁ω₃⁻¹ω₂ is the branch** - the other
#    order, ω₂ω₃⁻¹ω₁, is measured to be exactly −g (they are transposes of each other up
#    to sign, so this is the orientation of the quaternionic triple, not a free choice).
#    Symmetry and positive-definiteness of the result are **not** imposed; they are
#    oracles. Pinned against `semiflat_metric` at large R|Z| on SU(2).
# 8. Ooguri–Vafa: the Gibbons–Hawking potential is V = (R/4π)[−log(R²|z|²/μ²)] +
#    (R/π) Σ_{n≥1} cos(nθ_e) K₀(2πRn|z|), i.e. the additive constant is carried by the
#    cutoff μ. We take μ = 1 in the units of Z; the constant is invisible in the twistor
#    comparison, which is done on d log 𝒳 (oracle: `ooguri_vafa_xi` vs `solve_gmn`).

using QuadGK: quadgk

# ── the torus data ────────────────────────────────────────────────────────────────────

"""
    GMNTorus{F}

The input data of the finite-radius Gaiotto–Moore–Neitzke equations: the BPS charges
`charges[a]` (integer vectors in some lattice basis) with DT invariants `omega[a]`,
central charges `Z[a]` and quadratic refinement `sigma[a]`, the Gram matrix `pairing`
of the basis, the compactification radius `R`, and the torus angles `theta` (one per
basis charge, so ``θ_γ = Σ_i γ_i θ_i``).

Build with [`gmn_torus`](@ref); solve with [`solve_gmn`](@ref).
"""
struct GMNTorus{F}
    charges::Vector{Vector{Int}}
    omega::Vector{Int}
    Z::Vector{Complex{F}}
    sigma::Vector{Int}
    pairing::Matrix{Int}
    R::F
    theta::Vector{F}
    basisZ::Vector{Complex{F}}
    rays::Vector{Complex{F}}
    state_pairing::Matrix{Int}
    angles::Vector{F}
end

n_states(t::GMNTorus) = length(t.Z)

"""
    n_charges(t::GMNTorus) -> Int

The rank of the charge lattice.
"""
n_charges(t::GMNTorus) = size(t.pairing, 1)

"""
    radius(t::GMNTorus) -> Real

The compactification radius ``R`` of the circle.
"""
radius(t::GMNTorus) = t.R

"""
    torus_angles(t::GMNTorus) -> Vector

The torus angles ``θ_i``, one per basis charge.
"""
torus_angles(t::GMNTorus) = t.theta

central_charges(t::GMNTorus) = t.Z

"""
    su2_pairing() -> Matrix{Int}

The physical Dirac pairing ``⟨(n_m,n_e),(n_m',n_e')⟩ = n_m n_e' − n_e n_m'`` as a Gram
matrix, i.e. the pairing of the charges [`su2_bps_states`](@ref) returns. Pass it as the
`pairing` argument of [`gmn_torus`](@ref).
"""
su2_pairing() = [0 1; -1 0]

"""
    gmn_torus(states::AbstractVector{<:BPSState}, pairing::AbstractMatrix{<:Integer};
              R = 1, theta = nothing, sigma = -1, basis_Z = nothing) -> GMNTorus
    gmn_torus(sp::BPSSpectrum; kwargs...) -> GMNTorus

Assemble the finite-radius GMN data from a BPS spectrum. `pairing` is the Gram matrix
of the lattice basis the state charges are written in ([`su2_pairing`](@ref) for
[`su2_bps_states`](@ref), [`signed_pairing`](@ref) for [`bps_spectrum`](@ref) - the
latter is filled in automatically by the `BPSSpectrum` method). `R` is the
compactification radius, `theta` the torus angles (default all zero, one per basis
charge) and `sigma` the quadratic refinement, either one sign for every state or one
per state (default `-1`, the `log(1 + 𝒳)` convention of the conformal-limit TBA layer).

`basis_Z` gives the central charges of the *lattice basis*, which the metric needs even
for directions carrying no BPS state (the Ooguri-Vafa case: rank 2, one state). By
default they are recovered from the states, which requires the state charges to span
the lattice.
"""
function gmn_torus(states::AbstractVector{BPSState{F}},
                   pairing::AbstractMatrix{<:Integer};
                   R = 1, theta = nothing, sigma = -1, basis_Z = nothing) where {F}
    n = length(states)
    n ≥ 1 || throw(TBAError("no BPS states - there is no GMN system to solve"))
    r = size(pairing, 1)
    size(pairing, 2) == r || throw(Resurgence.InvalidArgument(
        "the pairing must be a square Gram matrix, got size $(size(pairing))"))
    cs = [charge(s) for s in states]
    all(length(c) == r for c in cs) || throw(Resurgence.InvalidArgument(
        "charge length does not match the rank $r of the pairing matrix"))
    Rf = F(R)
    Rf > 0 || throw(Resurgence.InvalidArgument("the radius R must be positive, got $R"))
    th = theta === nothing ? zeros(F, r) : F[F(x) for x in theta]
    length(th) == r || throw(Resurgence.InvalidArgument(
        "expected $r torus angles, got $(length(th))"))
    sg = sigma isa Integer ? fill(Int(sigma), n) : Int[Int(x) for x in sigma]
    length(sg) == n || throw(Resurgence.InvalidArgument(
        "expected $n quadratic-refinement signs, got $(length(sg))"))
    all(s -> s == 1 || s == -1, sg) || throw(Resurgence.InvalidArgument(
        "the quadratic refinement must be ±1, got $sg"))
    Z = Complex{F}[s.central_charge for s in states]
    any(iszero, Z) && throw(TBAError("a BPS state has vanishing central charge - its " *
                                     "ray is undefined; move off the singular locus"))
    P = Matrix{Int}(pairing)
    sp = [transpose(cs[a]) * P * cs[b] for a in 1:n, b in 1:n]
    bZ = basis_Z === nothing ? _infer_basis_Z(cs, Z, r) :
         Complex{F}[Complex{F}(z) for z in basis_Z]
    length(bZ) == r || throw(Resurgence.InvalidArgument(
        "expected $r basis central charges, got $(length(bZ))"))
    for a in 1:n
        isapprox(sum(cs[a][i] * bZ[i] for i in 1:r), Z[a];
                 rtol = sqrt(eps(F)), atol = sqrt(eps(F))) ||
            throw(Resurgence.InvalidArgument(
                "basis_Z is inconsistent with the spectrum: state $a has charge " *
                "$(cs[a]) and Z = $(Z[a]), but the basis gives " *
                "$(sum(cs[a][i] * bZ[i] for i in 1:r))"))
    end
    GMNTorus{F}(cs, [omega(s) for s in states], Z, sg, P, Rf, th, bZ,
                [-z / abs(z) for z in Z], sp,
                F[sum(c[i] * th[i] for i in 1:r) for c in cs])
end

# basis central charges from the states: solve M zs = Z with M the charge matrix. Only
# possible when the charges span the lattice - otherwise the metric would silently be
# built on an arbitrary Z in the missing directions, so we refuse and ask for basis_Z.
function _infer_basis_Z(cs::Vector{Vector{Int}}, Z::Vector{Complex{F}}, r::Int) where {F}
    M = Matrix{F}(reduce(hcat, cs)')
    LinearAlgebra.rank(M) == r || throw(Resurgence.InvalidArgument(
        "the BPS charges span a rank-$(LinearAlgebra.rank(M)) sublattice of the " *
        "rank-$r lattice, so the basis central charges are not determined - pass " *
        "them as `basis_Z`"))
    Complex{F}.(M \ Z)
end

gmn_torus(sp::BPSSpectrum; kwargs...) =
    gmn_torus(sp.states, signed_pairing(sp.basis); kwargs...)

"""
    GMNSolution{F}

A converged solution of the finite-radius GMN equations: the `torus`, the shared
rapidity `grid` (``ζ = e^{-θ} e_γ`` on each BPS ray), the sampled values
`gplus[k,a] =` ``\\log 𝒳_{γ_a}`` on its own ray and `gminus[k,a]` on the opposite ray,
the iteration count and the final sup-norm `residual`. Evaluate anywhere off the grid
with [`xi_value`](@ref).
"""
struct GMNSolution{F}
    torus::GMNTorus{F}
    grid::Vector{F}
    gplus::Matrix{Complex{F}}
    gminus::Matrix{Complex{F}}
    iterations::Int
    residual::F
end

n_iterations(sol::GMNSolution) = sol.iterations
residual(sol::GMNSolution) = sol.residual
n_states(sol::GMNSolution) = n_states(sol.torus)

# ── the solver ────────────────────────────────────────────────────────────────────────

# log(1 − σ𝒳) on the source rays; for σ = −1 these are tba.jl's _lplus / _lminus
_gmn_lplus(g, σ) = log1p(-σ * exp(g))
_gmn_lminus(g, σ) = log1p(-σ * exp(-g))

# the semiflat exponent −2πR|Z| cosh θ + iθ_γ, and its ζ → −ζ partner
_semiflat_plus(R::F, m, θ, ang) where {F} = -2 * F(π) * R * m * cosh(θ) + im * ang
_semiflat_minus(R::F, m, θ, ang) where {F} = 2 * F(π) * R * m * cosh(θ) + im * ang

# symmetric window: |𝒳^sf| = exp(−2πR m cosh θ) drops below `tol` outside it
function _gmn_window(R::F, mmin::F, tol::F) where {F}
    x = -log(tol) / (2 * F(π) * R * mmin)
    θmax = x > 1 ? acosh(x) : one(F)
    -θmax, θmax
end

"""
    solve_gmn(t::GMNTorus; window = nothing, n_points = nothing, tol = nothing,
              maxiter = 500, relax = 1, seed = nothing) -> GMNSolution

Solve the finite-radius Gaiotto–Moore–Neitzke equations of `t` by fixed-point iteration
from the semiflat seed ``\\log 𝒳_γ = πR Z_γ/ζ + iθ_γ + πR \\bar Z_γ ζ``. `window =
(θmin, θmax)` is the shared rapidity window (default symmetric, set by where the
semiflat exponential drops below `tol`), `n_points` the grid size, `relax ∈ (0,1]`
under-relaxes the update, and `seed` warm-starts from an earlier solution on an
identical grid. Throws [`TBAError`](@ref) if the sup-norm change has not fallen below
`tol` within `maxiter` sweeps.

Precision follows the element type of the central charges.
"""
function solve_gmn(t::GMNTorus{F}; window = nothing, n_points = nothing,
                   tol = nothing, maxiter::Integer = 500, relax::Real = 1,
                   seed::Union{Nothing,GMNSolution} = nothing) where {F}
    0 < relax ≤ 1 || throw(Resurgence.InvalidArgument(
        "relax must be in (0, 1], got $relax"))
    n = n_states(t)
    m = abs.(t.Z)
    tol = tol === nothing ? sqrt(eps(F)) * F(1 // 100) : F(tol)
    θmin, θmax = window === nothing ? _gmn_window(t.R, minimum(m), tol) :
                 (F(window[1]), F(window[2]))
    θmin < θmax || throw(TBAError("empty rapidity window [$θmin, $θmax]"))
    np = n_points === nothing ? max(101, ceil(Int, 12 * (θmax - θmin)) + 1) : Int(n_points)
    grid = collect(range(θmin, θmax; length = np))
    w = _trapezoid(grid)

    # Kernel blocks. On a uniform grid A(t_k e_a, ±t_l e_b) = (ρ(±e_b) + e_a)/(ρ(±e_b) −
    # e_a) with ρ = t_l/t_k = e^{(k−l)dθ}, so each block is **Toeplitz**: 2np−1 numbers
    # instead of np². That is what makes the nine solves of a metric point - and the
    # ten-state weak-coupling SU(2) chamber - affordable. The same two vectors serve the
    # opposite ray as well, since e_a → −e_a swaps A and B.
    dθ = grid[2] - grid[1]
    ρs = [exp(d * dθ) for d in -(np - 1):(np - 1)]
    KA = Dict{Tuple{Int,Int},Vector{Complex{F}}}()
    KB = Dict{Tuple{Int,Int},Vector{Complex{F}}}()
    for a in 1:n, b in 1:n
        t.state_pairing[a, b] == 0 && continue
        ea, eb = t.rays[a], t.rays[b]
        KA[(a, b)] = Complex{F}[(ρ * eb + ea) / (ρ * eb - ea) for ρ in ρs]
        KB[(a, b)] = Complex{F}[(-ρ * eb + ea) / (-ρ * eb - ea) for ρ in ρs]
    end

    sf_p = Complex{F}[_semiflat_plus(t.R, m[a], θ, t.angles[a]) for θ in grid, a in 1:n]
    sf_m = Complex{F}[_semiflat_minus(t.R, m[a], θ, t.angles[a]) for θ in grid, a in 1:n]
    if seed === nothing
        gp = copy(sf_p)
        gm = copy(sf_m)
    else
        size(seed.gplus) == (np, n) || throw(Resurgence.InvalidArgument(
            "the seed solution was solved on a different grid - pass matching " *
            "`window` and `n_points`"))
        gp = copy(seed.gplus)
        gm = copy(seed.gminus)
    end
    newgp = similar(gp)
    newgm = similar(gm)
    Lp = similar(gp)
    Lm = similar(gp)

    res = F(Inf)
    iters = 0
    for it in 1:maxiter
        iters = it
        for b in 1:n
            σ = t.sigma[b]
            @views @. Lp[:, b] = w * _gmn_lplus(gp[:, b], σ)      # weights folded in
            @views @. Lm[:, b] = w * _gmn_lminus(gm[:, b], σ)
        end
        copyto!(newgp, sf_p)
        copyto!(newgm, sf_m)
        for a in 1:n, b in 1:n
            c = t.state_pairing[a, b]
            c == 0 && continue
            coef = t.omega[b] * c / (4 * F(π) * im)
            A, B = KA[(a, b)], KB[(a, b)]
            for k in 1:np
                accp = zero(Complex{F})
                accm = zero(Complex{F})
                @inbounds for l in 1:np
                    d = k - l + np
                    accp += A[d] * Lp[l, b] - B[d] * Lm[l, b]
                    accm += B[d] * Lp[l, b] - A[d] * Lm[l, b]
                end
                newgp[k, a] += coef * accp
                newgm[k, a] += coef * accm
            end
        end
        res = max(maximum(abs, newgp - gp), maximum(abs, newgm - gm))
        @. gp = gp + relax * (newgp - gp)
        @. gm = gm + relax * (newgm - gm)
        res < tol && return GMNSolution{F}(t, grid, gp, gm, iters, res)
    end
    throw(TBAError("the GMN iteration did not converge: sup-norm change $res after " *
                   "$maxiter sweeps (tol $tol) - pass a larger maxiter, or relax < 1"))
end

# ── evaluation off the grid ───────────────────────────────────────────────────────────

# log 𝒳_γ(ζ) for an arbitrary lattice charge γ, by one pass of the integral equation over
# the converged ray data. `s = ±1` picks a lateral boundary value on a source ray (the
# kernel pole contributes ±½ of the residue), `s = 0` the principal value.
function _log_xi(sol::GMNSolution{F}, γ::AbstractVector{<:Integer}, ζ::Number,
                 s::Integer) where {F}
    t = sol.torus
    length(γ) == n_charges(t) || throw(Resurgence.InvalidArgument(
        "charge $γ does not have the lattice rank $(n_charges(t))"))
    iszero(ζ) && throw(TBAError("𝒳_γ has an essential singularity at ζ = 0"))
    n = n_states(t)
    grid, gp, gm = sol.grid, sol.gplus, sol.gminus
    ts = exp.(-grid)
    w = _trapezoid(grid)
    dθ = grid[2] - grid[1]
    Zγ = _charge_Z(t, γ)
    angγ = sum(γ[i] * t.theta[i] for i in 1:n_charges(t))
    acc = Complex{F}(F(π) * t.R * (Zγ / ζ + conj(Zγ) * ζ) + im * angγ)
    tabs = abs(ζ)
    θt = -log(tabs)
    for b in 1:n
        c = transpose(γ) * t.pairing * t.charges[b]
        c == 0 && continue
        σb = t.sigma[b]
        coef = t.omega[b] * c / (4 * F(π) * im)
        for ε in (1, -1)
            ej = ε * t.rays[b]
            col = ε == 1 ? view(gp, :, b) : view(gm, :, b)
            L = ε == 1 ? _gmn_lplus.(col, σb) : _gmn_lminus.(col, σb)
            onray = abs(angle(ζ / ej)) < sqrt(eps(F))
            if !onray
                accε = zero(Complex{F})
                @inbounds for k in eachindex(ts)
                    accε += w[k] * _tba_kernel(ζ, ts[k] * ej) * L[k]
                end
                acc += ε * coef * accε
            else
                gt = _interp_ray(grid, col, θt)
                Lt = ε == 1 ? _gmn_lplus(gt, σb) : _gmn_lminus(gt, σb)
                accε = zero(Complex{F})
                @inbounds for k in eachindex(ts)
                    if abs(grid[k] - θt) < dθ / 2
                        gl = _interp_ray(grid, col, θt - dθ)
                        gr = _interp_ray(grid, col, θt + dθ)
                        Ll = ε == 1 ? _gmn_lplus(gl, σb) : _gmn_lminus(gl, σb)
                        Lr = ε == 1 ? _gmn_lplus(gr, σb) : _gmn_lminus(gr, σb)
                        accε += w[k] * (-2) * (Lr - Ll) / (2 * dθ)
                    else
                        accε += w[k] * _tba_kernel(ζ, ts[k] * ej) * (L[k] - Lt)
                    end
                end
                accε += Lt * _pv_kernel_integral(grid, tabs)
                acc += ε * coef * accε
                acc += s * t.omega[b] * c * ε * Lt / 2
            end
        end
    end
    acc
end

"""
    charge_central_charge(t::GMNTorus, γ::AbstractVector{<:Integer}) -> Complex

The central charge ``Z_γ = Σ_i γ_i Z_{e_i}`` of an arbitrary lattice charge, linear in
`γ` over the basis central charges.
"""
function charge_central_charge(t::GMNTorus, γ::AbstractVector{<:Integer})
    length(γ) == n_charges(t) || throw(Resurgence.InvalidArgument(
        "charge $γ does not have the lattice rank $(n_charges(t))"))
    _charge_Z(t, γ)
end

_charge_Z(t::GMNTorus, γ) = sum(γ[i] * t.basisZ[i] for i in 1:n_charges(t))

"""
    xi_value(sol::GMNSolution, γ::AbstractVector{<:Integer}, ζ; side = :median) -> Complex

The Darboux coordinate ``𝒳_γ(ζ)`` of the lattice charge `γ`, evaluated by one pass of
the integral equation over the converged ray data. If `ζ` lies on a BPS ray the kernel
pole is principal-valued; `side = :plus`/`:minus` picks the lateral boundary value
instead, and the two sides differ by the Kontsevich–Soibelman factor
``(1 - σ(γ')𝒳_{γ'})^{Ω(γ')⟨γ,γ'⟩}``.
"""
xi_value(sol::GMNSolution, γ::AbstractVector{<:Integer}, ζ::Number;
         side::Symbol = :median) = exp(_log_xi(sol, γ, ζ, _side_sign(side)))

"""
    xi_value(sol::GMNSolution, i::Integer, ζ; side = :median) -> Complex

``𝒳_{γ_i}(ζ)`` for the `i`-th BPS state of the spectrum.
"""
xi_value(sol::GMNSolution, i::Integer, ζ::Number; side::Symbol = :median) =
    xi_value(sol, sol.torus.charges[i], ζ; side)

"""
    semiflat_xi(t::GMNTorus, γ::AbstractVector{<:Integer}, ζ) -> Complex

The uncorrected (semiflat) Darboux coordinate
``𝒳^{sf}_γ(ζ) = \\exp(πR Z_γ/ζ + iθ_γ + πR \\bar Z_γ ζ)``, the seed of
[`solve_gmn`](@ref) and its ``R|Z| → ∞`` limit.
"""
function semiflat_xi(t::GMNTorus{F}, γ::AbstractVector{<:Integer}, ζ::Number) where {F}
    r = n_charges(t)
    length(γ) == r || throw(Resurgence.InvalidArgument(
        "charge $γ does not have the lattice rank $r"))
    exp(_semiflat_log(t, γ, ζ))
end

# the semiflat exponent itself: log 𝒳^sf, which unlike semiflat_xi neither overflows at
# large |ζ| nor wraps the branch of the logarithm
function _semiflat_log(t::GMNTorus{F}, γ::AbstractVector{<:Integer}, ζ::Number) where {F}
    Zγ = _charge_Z(t, γ)
    angγ = sum(γ[i] * t.theta[i] for i in 1:n_charges(t))
    F(π) * t.R * (Zγ / ζ + conj(Zγ) * ζ) + im * angγ
end

"""
    instanton_correction(sol::GMNSolution, γ::AbstractVector{<:Integer}, ζ) -> Complex

The BPS-instanton correction ``\\log 𝒳_γ(ζ) - \\log 𝒳^{sf}_γ(ζ)``: everything the
integral equation adds to the semiflat answer. Computed on the logarithms, so it stays
meaningful at large ``|ζ|`` where ``𝒳^{sf}`` itself overflows.
"""
instanton_correction(sol::GMNSolution, γ::AbstractVector{<:Integer}, ζ::Number) =
    _log_xi(sol, γ, ζ, 0) - _semiflat_log(sol.torus, γ, ζ)

# ── the metric ────────────────────────────────────────────────────────────────────────
#
# The moduli space of a rank-2 charge lattice is four-dimensional: the complex modulus u
# and the two torus angles. Its coordinates here are always ordered
# (Re u, Im u, θ₁, θ₂), θ_i being the angle of the i-th basis charge.

"""
    MetricPoint{F}

The nine GMN solutions a metric evaluation needs: the base point and the two shifted
solutions per real coordinate, all on a **shared** rapidity grid so the central
differences are not polluted by a moving grid. Build with [`metric_point`](@ref);
consume with [`holomorphic_symplectic_form`](@ref), [`symplectic_expansion`](@ref) and
[`hk_metric`](@ref).
"""
struct MetricPoint{F}
    base::GMNSolution{F}
    plus::Vector{GMNSolution{F}}
    minus::Vector{GMNSolution{F}}
    steps::Vector{F}
end

n_charges(mp::MetricPoint) = n_charges(mp.base.torus)
radius(mp::MetricPoint) = radius(mp.base.torus)

"""
    metric_point(f, u::Number, theta::AbstractVector; h = 1e-3, htheta = h,
                 window = nothing, n_points = nothing, kwargs...) -> MetricPoint

Solve the GMN equations at a moduli-space point and at its eight central-difference
neighbours. `f(u, theta) -> GMNTorus` builds the BPS data at a given modulus and set of
torus angles - for pure ``SU(2)`` that is `(u, θ) -> su2_torus(sw, u, θ; R)`, see
[`su2_torus`](@ref). `h` is the step in `u` (both real and imaginary directions) and
`htheta` the step in the angles. Remaining keywords go to [`solve_gmn`](@ref); the
shared `window`/`n_points` default to the base point's own choice, slightly widened.
"""
function metric_point(f, u::Number, theta::AbstractVector; h = 1e-3, htheta = h,
                      window = nothing, n_points = nothing, kwargs...)
    t0 = f(u, theta)
    r = n_charges(t0)
    r == 2 || throw(Resurgence.InvalidArgument(
        "the metric layer needs a rank-2 charge lattice, got rank $r"))
    length(theta) == 2 || throw(Resurgence.InvalidArgument(
        "expected 2 torus angles, got $(length(theta))"))
    F = typeof(radius(t0))
    tol = get(kwargs, :tol, sqrt(eps(F)) * F(1 // 100))
    if window === nothing
        θmin, θmax = _gmn_window(t0.R, minimum(abs, t0.Z), F(tol))
        window = (θmin - one(F), θmax + one(F))
    end
    np = n_points === nothing ?
         max(101, ceil(Int, 12 * (window[2] - window[1])) + 1) : Int(n_points)
    base = solve_gmn(t0; window, n_points = np, kwargs...)
    hu, hθ = F(h), F(htheta)
    shifts = ((hu, 0), (im * hu, 0), (0, [hθ, zero(F)]), (0, [zero(F), hθ]))
    plus = GMNSolution{F}[]
    minus = GMNSolution{F}[]
    for (μ, (du, dθ)) in enumerate(shifts)
        for (store, s) in ((plus, 1), (minus, -1))
            up = μ ≤ 2 ? u + s * du : u
            θp = μ ≤ 2 ? collect(theta) : collect(theta) + s * dθ
            push!(store, solve_gmn(f(up, θp); window, n_points = np, seed = base,
                                   kwargs...))
        end
    end
    MetricPoint{F}(base, plus, minus, F[hu, hu, hθ, hθ])
end

# ∂_μ log 𝒳_{e_i}(ζ) by central differences: a 4 × 2 complex matrix
function _dlog_xi(mp::MetricPoint{F}, ζ::Number) where {F}
    D = Matrix{Complex{F}}(undef, 4, 2)
    for μ in 1:4, i in 1:2
        γ = [j == i ? 1 : 0 for j in 1:2]
        D[μ, i] = (_log_xi(mp.plus[μ], γ, ζ, 0) - _log_xi(mp.minus[μ], γ, ζ, 0)) /
                  (2 * mp.steps[μ])
    end
    D
end

"""
    holomorphic_symplectic_form(mp::MetricPoint, ζ) -> Matrix

The twistor-fibre holomorphic symplectic form
``ϖ(ζ) = (1/4π²R⟨γ₁,γ₂⟩) \\, d\\log 𝒳_{γ₁} ∧ d\\log 𝒳_{γ₂}`` as a 4×4 antisymmetric
matrix in the coordinates ``(\\mathrm{Re}\\,u, \\mathrm{Im}\\,u, θ₁, θ₂)``. The basis
must be unimodular symplectic (``⟨γ₁,γ₂⟩ = ±1``).
"""
function holomorphic_symplectic_form(mp::MetricPoint{F}, ζ::Number) where {F}
    p = mp.base.torus.pairing[1, 2]
    abs(p) == 1 || throw(Resurgence.InvalidArgument(
        "the lattice basis must be unimodular symplectic (⟨γ₁,γ₂⟩ = ±1), got $p"))
    D = _dlog_xi(mp, ζ)
    pref = 1 / (4 * F(π)^2 * radius(mp) * p)
    ϖ = Matrix{Complex{F}}(undef, 4, 4)
    for μ in 1:4, ν in 1:4
        ϖ[μ, ν] = pref * (D[μ, 1] * D[ν, 2] - D[ν, 1] * D[μ, 2])
    end
    ϖ
end

# sample points for the Laurent fit: phases kept clear of every BPS ray of the base
# point (where 𝒳 jumps), at two radii so the three powers separate cleanly
function _fit_zetas(mp::MetricPoint{F}) where {F}
    rays = mp.base.torus.rays
    forbidden = vcat([angle(e) for e in rays], [angle(-e) for e in rays])
    good = F[]
    for k in 0:47
        φ = F(2π) * k / 48
        d = minimum(abs(rem2pi(φ - f, RoundNearest)) for f in forbidden)
        d > F(0.2) && push!(good, φ)
        length(good) == 3 && break
    end
    length(good) == 3 || throw(TBAError(
        "could not find three ζ phases clear of the BPS rays - the spectrum's rays " *
        "fill the plane too densely for the Laurent fit"))
    Complex{F}[F(ρ) * cis(φ) for ρ in (F(0.8), F(1.25)) for φ in good]
end

"""
    symplectic_expansion(mp::MetricPoint; zetas = nothing)
        -> (; omega_plus, omega_3, omega_minus, residual)

Fit ``ϖ(ζ) = (-i/2ζ) ω_+ + ω_3 - (i/2) ζ ω_-`` from samples of
[`holomorphic_symplectic_form`](@ref). That only three powers of ``ζ`` occur is the
twistor theorem of [GMN10] - the `residual` of the least-squares fit is therefore a
numerical test that the ``𝒳_γ`` really are Darboux coordinates, and the returned
``ω_i`` must come out real.
"""
function symplectic_expansion(mp::MetricPoint{F}; zetas = nothing) where {F}
    ζs = zetas === nothing ? _fit_zetas(mp) : Complex{F}[Complex{F}(z) for z in zetas]
    length(ζs) ≥ 3 || throw(Resurgence.InvalidArgument(
        "the Laurent fit needs at least 3 ζ samples, got $(length(ζs))"))
    V = hcat(1 ./ ζs, ones(Complex{F}, length(ζs)), ζs)
    B = reduce(vcat, [transpose(vec(holomorphic_symplectic_form(mp, ζ))) for ζ in ζs])
    C = V \ B
    scale = max(maximum(abs, B), eps(F))
    residual = maximum(abs, V * C - B) / scale
    ωp = 2 * im * reshape(C[1, :], 4, 4)
    ω3 = reshape(C[2, :], 4, 4)
    ωm = 2 * im * reshape(C[3, :], 4, 4)
    (omega_plus = ωp, omega_3 = ω3, omega_minus = ωm, residual = residual)
end

# the hyperkähler triple, then g = ω₁ω₃⁻¹ω₂ (ledger item 7)
function _hk_metric_raw(mp::MetricPoint; zetas = nothing)
    e = symplectic_expansion(mp; zetas)
    ω1 = (e.omega_plus + e.omega_minus) / 2
    ω2 = (e.omega_plus - e.omega_minus) / (2 * im)
    g = ω1 * (e.omega_3 \ ω2)
    (g = g, omega_1 = ω1, omega_2 = ω2, omega_3 = e.omega_3, residual = e.residual)
end

"""
    hk_metric(mp::MetricPoint; zetas = nothing) -> Matrix
    hk_metric(f, u, theta; kwargs...) -> Matrix

The instanton-corrected hyperkähler metric on the moduli space, as a real symmetric
4×4 matrix in the coordinates ``(\\mathrm{Re}\\,u, \\mathrm{Im}\\,u, θ₁, θ₂)``. Built
from the hyperkähler triple of [`symplectic_expansion`](@ref) via
``g = ω_1 ω_3^{-1} ω_2``, which follows from ``J_i = g^{-1}ω_i`` and the quaternionic
relation among the triple.

That the result is real, symmetric and positive definite is a consequence of the
construction, not something imposed here - see [`hk_diagnostics`](@ref).
"""
function hk_metric(mp::MetricPoint; zetas = nothing)
    g = _hk_metric_raw(mp; zetas).g
    real.((g + transpose(g)) / 2)
end

hk_metric(f, u::Number, theta::AbstractVector; zetas = nothing, kwargs...) =
    hk_metric(metric_point(f, u, theta; kwargs...); zetas)

"""
    hk_diagnostics(mp::MetricPoint; zetas = nothing)
        -> (; residual, imaginary, asymmetry, eigenvalues)

Quality measures of a metric evaluation, all of which are oracles rather than knobs:
the Laurent-fit `residual` of [`symplectic_expansion`](@ref), the largest relative
`imaginary` part of the hyperkähler triple, the relative `asymmetry` of the raw
``ω_1ω_3^{-1}ω_2``, and the `eigenvalues` of the symmetrized metric (all positive for a
Riemannian metric).
"""
function hk_diagnostics(mp::MetricPoint; zetas = nothing)
    raw = _hk_metric_raw(mp; zetas)
    ims = maximum(maximum(abs, imag.(ω)) / max(maximum(abs, real.(ω)), eps())
                  for ω in (raw.omega_1, raw.omega_2, raw.omega_3))
    g = raw.g
    asym = maximum(abs, g - transpose(g)) / max(maximum(abs, g), eps())
    (residual = raw.residual, imaginary = ims, asymmetry = asym,
     eigenvalues = sort(real.(LinearAlgebra.eigvals(real.((g + transpose(g)) / 2)))))
end

"""
    semiflat_metric(t::GMNTorus, dZdu::AbstractVector) -> Matrix

The semiflat (uncorrected) metric

```math
g^{sf} = R\\,\\mathrm{Im}\\,τ\\,|da|^2
       + \\frac{1}{4π^2R}(\\mathrm{Im}\\,τ)^{-1}|dθ_1 - τ\\,dθ_2|^2,
```

with ``a = Z_{γ_2}``, ``a_D = Z_{γ_1}`` and ``τ = a_D'(u)/a'(u)`` built from the two
basis-period derivatives `dZdu`. This is the ``R|Z| → ∞`` limit of [`hk_metric`](@ref):
the BPS instanton corrections it drops are ``O(e^{-2πR|Z|})``.

For pure ``SU(2)`` the derivatives come from [`sw_periods`](@ref) - see
[`su2_torus`](@ref).
"""
function semiflat_metric(t::GMNTorus{F}, dZdu::AbstractVector) where {F}
    n_charges(t) == 2 || throw(Resurgence.InvalidArgument(
        "the semiflat metric needs a rank-2 charge lattice"))
    length(dZdu) == 2 || throw(Resurgence.InvalidArgument(
        "expected 2 period derivatives (dZ_{γ₁}/du, dZ_{γ₂}/du), got $(length(dZdu))"))
    aD′, a′ = Complex{F}(dZdu[1]), Complex{F}(dZdu[2])
    iszero(a′) && throw(TBAError("da/du vanishes - the special coordinate degenerates"))
    τ = aD′ / a′
    imτ = imag(τ)
    imτ > 0 || throw(TBAError("Im τ = $imτ ≤ 0 - the semiflat metric is not " *
                              "positive definite; check the period orientation"))
    R = t.R
    gu = R * imτ * abs2(a′)
    c = 1 / (4 * F(π)^2 * R)
    g = zeros(F, 4, 4)
    g[1, 1] = g[2, 2] = gu
    g[3, 3] = c / imτ
    g[3, 4] = g[4, 3] = -c * real(τ) / imτ
    g[4, 4] = c * (real(τ)^2 / imτ + imτ)
    g
end

"""
    su2_torus(sw::SeibergWittenSU2, u::Number, theta::AbstractVector; R = 1,
              chamber = :auto, tower = 4, sigma = -1) -> GMNTorus

The finite-radius GMN data of pure ``SU(2)`` at Coulomb modulus `u`: the BPS states of
[`su2_bps_states`](@ref) in the physical charge basis, the Dirac pairing
[`su2_pairing`](@ref), and the basis central charges ``(a_D, a)`` from
[`sw_periods`](@ref) - which are needed explicitly because the electric direction
``(0,1)`` carries no BPS state of its own.

Pass `(u, θ) -> su2_torus(sw, u, θ; R)` to [`metric_point`](@ref) to get the metric on
the SU(2) moduli space.
"""
function su2_torus(sw::SeibergWittenSU2, u::Number, theta::AbstractVector; R = 1,
                   chamber::Symbol = :auto, tower::Integer = 4, sigma = -1)
    states = su2_bps_states(sw, u; chamber, tower)
    a, aD = sw_periods(sw, u)
    gmn_torus(states, su2_pairing(); R, theta, sigma, basis_Z = [aD, a])
end

# ── Ooguri–Vafa: the one-hypermultiplet case, solvable in closed form ─────────────────
#
# With a single BPS state γ_s the equation degenerates: ⟨γ_s,γ_s⟩ = 0, so 𝒳_{γ_s} is
# exactly semiflat and the equation for every other charge is a *single explicit
# quadrature* of a known function - no fixed point at all. That gives two independent
# checks of the solver: adaptive quadrature of the same integral, and, in the ζ → ∞
# limit, a closed-form Bessel-K instanton sum. This is the Ooguri–Vafa geometry [OV96].

# K₀(x) = ∫₀^∞ e^{−x cosh t} dt, in-file (as the Carlson forms are in sw_curve.jl) so the
# package needs no special-function dependency
function _besselk0(x::F) where {F<:Real}
    x > 0 || throw(Resurgence.InvalidArgument("K₀ needs a positive argument, got $x"))
    val, _ = quadgk(t -> exp(-x * cosh(t)), zero(F), F(Inf); rtol = sqrt(eps(F)))
    val
end

# the single BPS state of an Ooguri-Vafa torus, or an error
function _ov_state(t::GMNTorus)
    n_states(t) == 1 || throw(Resurgence.InvalidArgument(
        "the Ooguri-Vafa case needs exactly one BPS state, got $(n_states(t))"))
    n_charges(t) == 2 || throw(Resurgence.InvalidArgument(
        "the Ooguri-Vafa case needs a rank-2 charge lattice"))
    1
end

"""
    ooguri_vafa_xi(t::GMNTorus, γ::AbstractVector{<:Integer}, ζ) -> Complex

``𝒳_γ(ζ)`` for a torus carrying a **single** BPS state (the Ooguri-Vafa geometry),
computed by adaptive quadrature of the GMN integral rather than by iteration - possible
because the lone state's own coordinate is exactly semiflat, so the equation for every
other charge is one explicit integral. Independent of [`solve_gmn`](@ref) in grid,
weights, window and iteration, and therefore its oracle. `ζ` must lie off the two BPS
rays.
"""
function ooguri_vafa_xi(t::GMNTorus{F}, γ::AbstractVector{<:Integer}, ζ::Number) where {F}
    s = _ov_state(t)
    c = transpose(collect(γ)) * t.pairing * t.charges[s]
    base = semiflat_xi(t, γ, ζ)
    c == 0 && return base
    e = t.rays[s]
    for sgn in (1, -1)
        abs(angle(ζ / (sgn * e))) > sqrt(eps(F)) || throw(TBAError(
            "ζ = $ζ lies on a BPS ray, where the quadrature is only defined as a " *
            "principal value - move off the ray"))
    end
    b = 2 * F(π) * t.R * abs(t.Z[s])
    ang = t.angles[s]
    σ = t.sigma[s]
    integrand(θ) = _tba_kernel(ζ, exp(-θ) * e) * log1p(-σ * exp(-b * cosh(θ) + im * ang)) -
                   _tba_kernel(ζ, -exp(-θ) * e) * log1p(-σ * exp(-b * cosh(θ) - im * ang))
    # the integrand decays like exp(−b cosh θ): integrate over the window where it is
    # above machine precision rather than over the whole line, where the kernel's
    # exponential growth would make the endpoint transform overflow
    x = -log(eps(F)) / b
    θmax = (x > 1 ? acosh(x) : one(F)) + one(F)
    val, _ = quadgk(integrand, -θmax, zero(F), θmax; rtol = sqrt(eps(F)))
    base * exp(t.omega[s] * c / (4 * F(π) * im) * val)
end

"""
    ooguri_vafa_instantons(t::GMNTorus, γ::AbstractVector{<:Integer}; nmax = 64) -> Real

The ``ζ → ∞`` limit of the instanton correction to ``\\log 𝒳_γ`` in the Ooguri-Vafa
geometry, in closed form:

```math
-\\frac{Ω⟨γ,γ_s⟩}{π} \\sum_{n≥1} \\frac{(-1)^{n+1}}{n}
    \\sin(n θ_{γ_s})\\, K_0(2πnR|Z_{γ_s}|).
```

Both kernels tend to ``-1`` as ``ζ → ∞``, leaving ``∫dθ\\, e^{-2πnR|Z|\\cosh θ} =
2K_0(2πnR|Z|)`` per instanton number. This is the Ooguri-Vafa instanton sum, and
comparing it with [`solve_gmn`](@ref) at large ``|ζ|`` is what pins the ``2πR`` of the
semiflat exponent (ledger item 2). The sum is truncated where the Bessel terms
underflow, at most `nmax` terms.
"""
function ooguri_vafa_instantons(t::GMNTorus{F}, γ::AbstractVector{<:Integer};
                                nmax::Integer = 64) where {F}
    s = _ov_state(t)
    c = transpose(collect(γ)) * t.pairing * t.charges[s]
    c == 0 && return zero(F)
    b = 2 * F(π) * t.R * abs(t.Z[s])
    ang = t.angles[s]
    acc = zero(F)
    for n in 1:nmax
        k = _besselk0(n * b)
        term = (isodd(n) ? one(F) : -one(F)) / n * sin(n * ang) * k
        acc += term
        k < eps(F) && break
    end
    -t.omega[s] * c * acc / F(π)
end

"""
    ooguri_vafa_torus(z::Number, theta::AbstractVector; R = 1, sigma = -1) -> GMNTorus

The Ooguri-Vafa local model: a rank-2 lattice with a single hypermultiplet of charge
``γ_s = (0,1)``, ``Ω = 1`` and central charge ``Z_{γ_s} = z``, whose dual basis charge
carries the one-loop period ``Z_{γ_1} = (z\\log z - z)/2πi``. The state becomes massless
at `z = 0`, where the semiflat metric degenerates and the instanton corrections resolve
it.
"""
function ooguri_vafa_torus(z::Number, theta::AbstractVector; R = 1, sigma = -1)
    iszero(z) && throw(Resurgence.InvalidArgument(
        "the Ooguri-Vafa modulus z = 0 is the singular point itself"))
    F = typeof(float(real(z)))
    Zs = Complex{F}(z)
    Zd = (Zs * log(Zs) - Zs) / (2 * F(π) * im)
    state = BPSState{F}([0, 1], Zs, 1)
    gmn_torus([state], su2_pairing(); R, theta, sigma, basis_Z = [Zd, Zs])
end

"""
    ooguri_vafa_period_derivatives(z::Number) -> Vector

The two basis-period `z`-derivatives of [`ooguri_vafa_torus`](@ref),
``(dZ_{γ_1}/dz, dZ_{γ_2}/dz) = (\\log z/2πi,\\ 1)``, ready for
[`semiflat_metric`](@ref). Their ratio ``τ = \\log z/2πi`` has
``\\mathrm{Im}\\,τ = -\\log|z|/2π > 0`` for ``|z| < 1``, which is the regime where the
local model is Riemannian.
"""
function ooguri_vafa_period_derivatives(z::Number)
    F = typeof(float(real(z)))
    [log(Complex{F}(z)) / (2 * F(π) * im), Complex{F}(1)]
end
