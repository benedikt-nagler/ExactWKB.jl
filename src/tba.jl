# TBA layer: exact quantum periods from the BPS spectrum alone. The conformal-limit
# Gaiotto–Moore–Neitzke integral equations (the form Ito–Mariño–Shu derived for
# polynomial oscillators) determine the Borel-summed Voros symbols non-perturbatively
# from nothing but the BPS data (Z_γ, Ω(γ), ⟨γ,γ'⟩) - no WKB series, no Padé. This is
# the independent cross-check of the whole series pipeline: wkb_expansion → voros →
# Borel–Padé lateral sums.
#
# ── Convention ledger (each item pinned by an oracle in test/test_tba.jl) ──────────
#
# 1. The equation. For every spectrum charge γ (and X_{−γ} = X_γ⁻¹ identically):
#      log X_γ(ħ) = Z_γ/ħ + Σ_{γ'} Ω(γ') ⟨γ,γ'⟩/(4πi) ∫_{ℓ_{γ'}} dħ'/ħ'
#                                        · (ħ'+ħ)/(ħ'−ħ) · log(1 + X_{γ'}(ħ')),
#    the sum running over ±(spectrum charges), with the BPS ray
#    ℓ_{γ'} = {ħ : Z_{γ'}/ħ ∈ ℝ_{<0}} - the ray where X_{γ'} decays (DDP ledger
#    item 4). The ±γ' pair is always integrated together on a shared grid: their
#    divergent kernel parts cancel exactly (the pair kernel decays at both ends),
#    which is what makes the conformal limit finite.
# 2. Identification: X_γ IS the summed Voros symbol V_γ of the DDP ledger, in the
#    same physical (signed-frame) charge basis: X evaluated off the rays equals
#    exp(v_{-1}/ħ + s(tail)) from voros_value. Pinned by the cubic and quartic
#    flagship oracles.
# 3. Boundary values on a ray. For ħ on a source ray the kernel pole gives
#      lateral(±) = PV ± Ω(γ')⟨γ,γ'⟩ · log(1 + X_{γ'}(ħ))/2,
#    so the plus/minus jump is the DDP connection formula
#    Δ log X_γ = Ω⟨γ,γ'⟩·log(1 + X_{γ'}) by construction; `side = :median` is the
#    principal value. The sign attaching :plus to +iπ is pinned by the jump oracle
#    against verify_ddp.
# 4. Rapidity parametrization: ħ = e^{−θ}·e_γ on ℓ_γ with e_γ = −Z_γ/|Z_γ|, so the
#    semiflat seed is log X_γ = −|Z_γ| e^{θ} (decaying, double-exponential UV falloff)
#    and the θ → −∞ end sits on the constant-Y-system plateau (for the cubic: the
#    golden ratio, Zamolodchikov's constant A₂ TBA).

"""
    TBASystem{F}

The data of the conformal-limit TBA equations extracted from a [`BPSSpectrum`](@ref):
per state its physical central charge `Z`, DT invariant `omega`, the state-level
intersection pairing `pairing[a,b] = ⟨γ_a, γ_b⟩`, and the BPS-ray direction
`rays[a] = −Z_a/|Z_a|`. Build with [`tba_system`](@ref); solve with
[`solve_tba`](@ref).
"""
struct TBASystem{F}
    spectrum::BPSSpectrum{F}
    Z::Vector{Complex{F}}
    omega::Vector{Int}
    pairing::Matrix{Int}
    rays::Vector{Complex{F}}
end

n_states(sys::TBASystem) = length(sys.Z)

"""
    tba_system(sp::BPSSpectrum) -> TBASystem

Extract the TBA data from a BPS spectrum: state central charges ``Z_γ = Σ c_j Z_j``
over the physical (signed-frame) basis charges, DT invariants ``Ω(γ)``, and the
state-level pairing ``⟨γ_a, γ_b⟩ = c_a^T P c_b`` with `P = signed_pairing(basis)`.
"""
function tba_system(sp::BPSSpectrum{F}) where {F}
    n_states(sp) ≥ 1 || throw(TBAError("the spectrum has no BPS states - " *
                                       "there is no TBA system to solve"))
    P = signed_pairing(sp.basis)
    cs = charges(sp)
    pairing = [transpose(cs[a]) * P * cs[b] for a in eachindex(cs), b in eachindex(cs)]
    Z = central_charges(sp)
    rays = [-z / abs(z) for z in Z]
    TBASystem{F}(sp, Z, [omega(s) for s in sp.states], pairing, rays)
end

"""
    TBASolution{F}

A converged solution of the TBA equations: the `system`, the shared rapidity `grid`
(uniform in θ, with ``ħ = e^{-θ} e_γ`` on each ray), the sampled values
`gplus[k,a] =` ``\\log X_{γ_a}`` on its own ray ``ℓ_{γ_a}`` and `gminus[k,a]` on the
opposite ray ``ℓ_{-γ_a}``, plus the iteration count and final sup-norm `residual`.
Evaluate anywhere off the grid with [`voros_value`](@ref).
"""
struct TBASolution{F}
    system::TBASystem{F}
    grid::Vector{F}
    gplus::Matrix{Complex{F}}
    gminus::Matrix{Complex{F}}
    iterations::Int
    residual::F
end

"""
    n_iterations(sol::TBASolution) -> Int

Number of fixed-point sweeps the TBA solver used.
"""
n_iterations(sol::TBASolution) = sol.iterations

"""
    residual(sol::TBASolution) -> Real

The final sup-norm change of `log X` between the last two sweeps.
"""
residual(sol::TBASolution) = sol.residual

# the pair kernel pieces: A(ħ, ħ') = (ħ'+ħ)/(ħ'−ħ), and the log-measure weight is
# carried by the uniform θ grid (dħ'/ħ' = dθ' after orientation)
_tba_kernel(ħ, ħp) = (ħp + ħ) / (ħp - ħ)

# trapezoid weights on the uniform grid
function _trapezoid(grid::Vector{F}) where {F}
    dθ = grid[2] - grid[1]
    w = fill(dθ, length(grid))
    w[1] = w[end] = dθ / 2
    w
end

# L⁺ = log(1 + X_γ) on ℓ_γ, L⁻ = log(1 + X_{−γ}) = log(1 + X_γ⁻¹) on ℓ_{−γ}
_lplus(g) = log1p(exp(g))
_lminus(g) = log1p(exp(-g))

# the folded ±γ_j pair contribution at target ħ (target charge index i), from the
# sampled sources: Ω_j ⟨γ_i,γ_j⟩/(4πi) Σ_k w_k [A(ħ, t_k e_j) Lp[k] − A(ħ, −t_k e_j) Lm[k]]
function _pair_sum(ħ, ts, w, ej, Lp::AbstractVector, Lm::AbstractVector)
    acc = zero(promote_type(typeof(ħ), eltype(Lp)))
    @inbounds for k in eachindex(ts)
        ħp = ts[k] * ej
        acc += w[k] * (_tba_kernel(ħ, ħp) * Lp[k] - _tba_kernel(ħ, -ħp) * Lm[k])
    end
    acc
end

function _default_window(mmin::F, tol::F) where {F}
    θmax = log(-F(6) / 5 * log(eps(F)) / mmin)
    θmin = log(tol) - 2
    θmin < θmax || throw(TBAError("empty rapidity window [$θmin, $θmax] - loosen " *
                                  "tol or pass window explicitly"))
    θmin, θmax
end

"""
    solve_tba(sp::BPSSpectrum; window = nothing, n_points = nothing, tol = 1e-10,
              maxiter = 500, relax = 1) -> TBASolution

Solve the conformal-limit TBA equations of the spectrum by fixed-point iteration
from the semiflat seed ``\\log X_γ = Z_γ/ħ``. `window = (θmin, θmax)` is the shared
rapidity window (default: UV end where the semiflat seed underflows, IR end set by
`tol`); `n_points` the grid size (default 12 points per unit of θ); `relax ∈ (0, 1]`
under-relaxes the update. Throws [`TBAError`](@ref) if the sup-norm change has not
fallen below `tol` within `maxiter` sweeps.

Precision follows the element type of the spectrum's central charges.
"""
function solve_tba(sp::BPSSpectrum{F}; window = nothing, n_points = nothing,
                   tol = nothing, maxiter::Integer = 500, relax::Real = 1) where {F}
    0 < relax ≤ 1 || throw(Resurgence.InvalidArgument(
        "relax must be in (0, 1], got $relax"))
    sys = tba_system(sp)
    n = n_states(sys)
    m = abs.(sys.Z)
    tol = tol === nothing ? sqrt(eps(F)) * F(1 // 100) : F(tol)
    θmin, θmax = window === nothing ? _default_window(minimum(m), tol) :
                 (F(window[1]), F(window[2]))
    np = n_points === nothing ? max(101, ceil(Int, 12 * (θmax - θmin)) + 1) :
         Int(n_points)
    grid = collect(range(θmin, θmax; length = np))
    w = _trapezoid(grid)
    ts = exp.(-grid)                               # t = e^{-θ}, ħ = ±t·e_j on ℓ_{±γ_j}

    # semiflat seed: on ℓ_{γ_a} log X = −m_a e^{θ}, on ℓ_{−γ_a} it is +m_a e^{θ}
    gp = Complex{F}[-m[a] * exp(θ) for θ in grid, a in 1:n]
    gm = Complex{F}[+m[a] * exp(θ) for θ in grid, a in 1:n]
    newgp = similar(gp)
    newgm = similar(gm)
    Lp = similar(gp)
    Lm = similar(gp)

    res = F(Inf)
    iters = 0
    for it in 1:maxiter
        iters = it
        @. Lp = _lplus(gp)
        @. Lm = _lminus(gm)
        for a in 1:n
            for k in eachindex(grid)
                ħ₊ = ts[k] * sys.rays[a]
                accp = sys.Z[a] / ħ₊
                accm = sys.Z[a] / (-ħ₊)
                for b in 1:n
                    b == a && continue
                    c = sys.pairing[a, b]
                    c == 0 && continue
                    coef = sys.omega[b] * c / (4 * F(π) * im)
                    accp += coef * _pair_sum(ħ₊, ts, w, sys.rays[b],
                                             view(Lp, :, b), view(Lm, :, b))
                    accm += coef * _pair_sum(-ħ₊, ts, w, sys.rays[b],
                                             view(Lp, :, b), view(Lm, :, b))
                end
                newgp[k, a] = accp
                newgm[k, a] = accm
            end
        end
        res = max(maximum(abs, newgp - gp), maximum(abs, newgm - gm))
        @. gp = gp + relax * (newgp - gp)
        @. gm = gm + relax * (newgm - gm)
        res < tol && return TBASolution{F}(sys, grid, gp, gm, iters, res)
    end
    throw(TBAError("the TBA iteration did not converge: sup-norm change $res after " *
                   "$maxiter sweeps (tol $tol) - pass a larger maxiter, or " *
                   "relax < 1"))
end

# linear interpolation of a sampled ray function at rapidity θ
function _interp_ray(grid::Vector{F}, col, θ) where {F}
    grid[1] ≤ θ ≤ grid[end] || throw(TBAError(
        "evaluation point at rapidity θ = $θ lies outside the solved window " *
        "[$(grid[1]), $(grid[end])] - re-solve with a wider window"))
    k = clamp(searchsortedlast(grid, θ), 1, length(grid) - 1)
    λ = (θ - grid[k]) / (grid[k + 1] - grid[k])
    (1 - λ) * col[k] + λ * col[k + 1]
end

_side_sign(side::Symbol) =
    side === :plus ? 1 :
    side === :minus ? -1 :
    side === :median ? 0 :
    throw(Resurgence.InvalidArgument(
        "side must be :plus, :minus or :median, got :$side"))

# closed-form ∫_{θmin}^{θmax} (t'+t)/(t'−t) dθ' (principal value across t' = t):
# with u = e^{−θ'} the measure dθ' = −du/u flips the limits, so the antiderivative
# 2·log|u−t| − log u is evaluated from u(θmax) up to u(θmin)
function _pv_kernel_integral(grid::Vector{F}, t) where {F}
    ua, ub = exp(-grid[1]), exp(-grid[end])
    (2 * log(abs(ua - t)) - log(ua)) - (2 * log(abs(ub - t)) - log(ub))
end

# log X_{γ_i}(ħ) evaluated anywhere inside the solved window. On a source BPS ray the
# kernel pole is handled by subtraction + the closed-form PV integral, and `s = ±1`
# adds the lateral ±iπ residue term (ledger item 3); s = 0 is the principal value.
function _log_x(sol::TBASolution{F}, i::Integer, ħ::Number, s::Integer) where {F}
    sys = sol.system
    n = n_states(sys)
    1 ≤ i ≤ n || throw(Resurgence.InvalidArgument(
        "state index $i is out of range 1:$n"))
    grid, gp, gm = sol.grid, sol.gplus, sol.gminus
    ts = exp.(-grid)
    w = _trapezoid(grid)
    dθ = grid[2] - grid[1]
    t = abs(ħ)
    θt = -log(t)
    acc = Complex{F}(sys.Z[i] / ħ)
    for b in 1:n
        b == i && continue
        c = sys.pairing[i, b]
        c == 0 && continue
        coef = sys.omega[b] * c / (4 * F(π) * im)
        for σ in (1, -1)
            ej = σ * sys.rays[b]
            L = σ == 1 ? _lplus.(view(gp, :, b)) : _lminus.(view(gm, :, b))
            # is ħ on this ray? (positive multiple of ej)
            onray = abs(angle(ħ / ej)) < sqrt(eps(F))
            if !onray
                accσ = zero(Complex{F})
                @inbounds for k in eachindex(ts)
                    accσ += w[k] * _tba_kernel(ħ, ts[k] * ej) * L[k]
                end
                acc += σ * coef * accσ
            else
                # pole at t' = t: subtract, integrate the kernel in closed form,
                # and add the lateral residue term per ledger item 3
                gt = _interp_ray(grid, view(σ == 1 ? gp : gm, :, b), θt)
                Lt = σ == 1 ? _lplus(gt) : _lminus(gt)
                accσ = zero(Complex{F})
                @inbounds for k in eachindex(ts)
                    if abs(grid[k] - θt) < dθ / 2
                        # the finite limit A·(L − Lt) → −2·dL/dθ at the pole
                        h = dθ
                        gl = _interp_ray(grid, view(σ == 1 ? gp : gm, :, b), θt - h)
                        gr = _interp_ray(grid, view(σ == 1 ? gp : gm, :, b), θt + h)
                        Ll = σ == 1 ? _lplus(gl) : _lminus(gl)
                        Lr = σ == 1 ? _lplus(gr) : _lminus(gr)
                        accσ += w[k] * (-2) * (Lr - Ll) / (2 * h)
                    else
                        accσ += w[k] * _tba_kernel(ħ, ts[k] * ej) * (L[k] - Lt)
                    end
                end
                accσ += Lt * _pv_kernel_integral(grid, t)
                acc += σ * coef * accσ
                acc += s * sys.omega[b] * c * σ * Lt / 2
            end
        end
    end
    acc
end

"""
    voros_value(sol::TBASolution, i::Integer, ħ; side = :median) -> Complex

The TBA-side value of the summed Voros symbol ``X_{γ_i}(ħ)`` of spectrum state `i`,
evaluated by one pass of the integral equation over the converged ray data. If `ħ`
lies on a BPS ray the kernel pole is principal-valued; `side = :plus`/`:minus` picks
the lateral boundary value instead, and the two sides differ by exactly the
Delabaere–Dillinger–Pham jump ``(1 + X_{γ'})^{Ω⟨γ,γ'⟩}``.

This is the non-perturbative counterpart of `voros_value(::VorosSymbol, ħ)`: the two
must agree wherever both converge - that identity is the flagship oracle of the TBA
layer.
"""
function voros_value(sol::TBASolution, i::Integer, ħ::Number; side::Symbol = :median)
    exp(_log_x(sol, i, ħ, _side_sign(side)))
end
