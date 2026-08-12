"""
    ExactWKB

Bridge package for the exact WKB / Stokes-graph layer of the cluster ecosystem.
Turns a one-dimensional Schrödinger problem ``ħ²ψ″ = Q(z)ψ`` (`Q = V − E`,
Iwaki–Nakanishi normalization) into the objects of exact WKB analysis: turning
points, the all-orders WKB (Riccati) recursion, quantum periods and Voros symbols,
and Stokes graphs with their saddle (BPS) data.

It depends on `Resurgence.jl` (the resurgence foundation - Borel/Padé/Laplace
summation of the Voros series) and on `ClusterAlgebras.jl` (the cluster core, live
since the DDP layer: Voros wall-crossing = y-mutation, see `src/ddp.jl`).
Resurgence names are **not** re-exported, and
AbstractAlgebra ring elements never leak into the public API - public types expose
plain numbers and `Resurgence.FormalSeries`.
"""
module ExactWKB

using PrecompileTools: @setup_workload, @compile_workload
import LinearAlgebra

# Foundation packages. Resurgence supplies FormalSeries + the summation toolbox;
# ClusterAlgebras supplies the y-mutation / seed machinery of the DDP layer.
import Resurgence
import ClusterAlgebras

export ExactWKBError, InvalidPotential, UnsupportedTurningPoint, ContourError,
       TracingFailed, NonGenericGraph, ChamberError, CoalescentTurningPoints,
       QuantizationError, TBAError, PeriodError
export SpectralCycle, spectral_cycles, uniform_cycles, quantum_period, kind, endpoints
export quantization_condition, spectral_determinant, wkb_eigenvalue,
       quantization_derivative, eigenvalue_sensitivity
export perturbative_b, instanton_a, verify_zjj, energy_splitting
export AbstractSchrodingerProblem, SchrodingerProblem, q_coefficients, energy, degree,
       variable, with_energy, q_derivative_at, q_taylor_at, asymptotic_directions
export RationalProblem, q_numerator, poles, pole_orders, n_finite_poles,
       mathieu_problem, puncture_indices, n_punctures,
       RingDomainWall, ring_domain_walls, pole_index
export TurningPoint, turning_points, simple_turning_points, location, order, is_simple
export WKBExpansion, wkb_expansion, s_odd_terms, evaluate_s_odd, even_odd_residual
export WKBDerivative, wkb_derivative
export period_integral, wkb_period, encircling_contour, period_derivative
export VorosSymbol, voros_symbol, classical_period, quantum_series, full_series
export WeberModel, weber_model, weber_problem, weber_index, merging_period,
       weber_voros_series, weber_voros_coefficient, weber_voros_sum, weber_log_gamma,
       weber_connection, weber_barrier_amplitude
export StokesLine, StokesGraph, stokes_graph, is_finite_line
export finite_lines, edges, n_infinite_lines, topology_signature,
       is_degenerate, turning_point_orders
export mass, endpoint
export Saddle, saddle_candidates, saddles, is_saddle, stokes_graph_family, central_charge
export voros_value, ddp_transform, verify_ddp, verify_ddp_mutation, ddp_seed,
       intersection_pairing
export IdealTriangulation, ideal_triangulation, triangulation_quiver, n_diagonals,
       diagonals, n_marked_points, n_boundaries, infinite_lines, boundary_lines,
       ray_exit_angles, flip, canonical_reorder, diagonal_core_paths,
       puncture_valence, selffolded_triangles, is_selffolded, selffold_arcs,
       TaggedArc, tagged_arcs
export PolygonDecomposition, polygon_decomposition, cells, cell_tp, cell_sizes,
       n_cells, is_triangulation, refinements, n_refinements
export ChargeBasis, charge_basis, charge_contour, bridge_seed, n_charges,
       central_charges, signs, physical_charges, signed_pairing
export signed_frame, verify_signed_frame, chamber_walls, reference_theta
export BPSState, BPSSpectrum, bps_spectrum, n_states, charges, charge, phase
export refined_skew, refined_dilog_word, refined_dt, refined_frame_map,
       refined_quiver_iso,
       verify_refined_wall_crossing, refined_ddp_jump, verify_refined_ddp
export TBASystem, TBASolution, tba_system, solve_tba, n_iterations, residual
export GMNTorus, GMNSolution, gmn_torus, solve_gmn, xi_value, semiflat_xi,
       su2_pairing, radius, torus_angles, charge_central_charge
export MetricPoint, metric_point, holomorphic_symplectic_form, symplectic_expansion,
       hk_metric, hk_diagnostics, semiflat_metric, su2_torus
export instanton_correction
export ooguri_vafa_xi, ooguri_vafa_instantons, ooguri_vafa_torus,
       ooguri_vafa_period_derivatives
export SeibergWittenSU2, sw_periods, sw_period_derivatives, quantum_sw_periods,
       sw_singularities,
       sw_monodromy, dynamical_scale, continue_periods, ms_wall, sw_chamber
export su2_bps_quiver, su2_bps_states, verify_su2_wall_crossing
export su2_refined_rays, verify_su2_refined_wall_crossing
export plot_stokes_graph, plot_triangulation

include("errors.jl")
include("potentials.jl")
include("rational_potentials.jl")
include("turning_points.jl")
include("wkb_recursion.jl")
include("periods.jl")
include("voros.jl")
include("weber.jl")
include("quantum_periods.jl")
include("stokes_graph.jl")
include("saddles.jl")
include("ddp.jl")
include("quantization.jl")
include("double_well.jl")
include("polygon_decomposition.jl")
include("triangulation.jl")
include("charge_lattice.jl")
include("signed_frame.jl")
include("bps.jl")
include("tba.jl")
include("refined.jl")
include("sw_curve.jl")
include("sw_bps.jl")
include("hyperkahler.jl")
include("show.jl")

"""
    plot_stokes_graph(g::StokesGraph; kwargs...)
    plot_stokes_graph(gs::AbstractVector{<:StokesGraph}; kwargs...)

Plot a Stokes graph (turning points, traced lines, highlighted saddles) or a θ-family
of them. Provided by the Makie extension - load a backend first (`using CairoMakie` or
`using GLMakie`). Without one this throws a hint.
"""
function plot_stokes_graph(args...; kwargs...)
    error("plot_stokes_graph requires a Makie backend - run `using CairoMakie` " *
          "(or `using GLMakie`) to load the ExactWKB Makie extension.")
end

"""
    plot_triangulation(g::StokesGraph, t::IdealTriangulation; kwargs...)

Plot the Stokes graph with its dual ideal triangulation overlaid: the boundary circle
with the marked points at the asymptotic directions, and the diagonals (labelled by
their turning-point pairs) as chords. Provided by the Makie extension - load a
backend first (`using CairoMakie` or `using GLMakie`).
"""
function plot_triangulation(args...; kwargs...)
    error("plot_triangulation requires a Makie backend - run `using CairoMakie` " *
          "(or `using GLMakie`) to load the ExactWKB Makie extension.")
end

# Precompile the flagship Stokes-graph chain on a cheap Float64 double well so the first
# interactive `stokes_graph`/`voros_symbol` call is warm.
@setup_workload begin
    @compile_workload begin
        prob = SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0])
        tps = turning_points(prob)
        w = wkb_expansion(prob; order = 2)
        voros_symbol(w, encircling_contour(tps[2], tps[3]))
        g = stokes_graph(prob; theta = 0.0)
        topology_signature(g)
        saddle_candidates(prob)
        # the cluster-bridge chain on the cheap cubic
        cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
        t = ideal_triangulation(stokes_graph(cubic; theta = 0.3))
        cb = charge_basis(cubic, t)
        bridge_seed(cb)
        signs(cb); physical_charges(cb); signed_pairing(cb)
        flip(t, 1)
        chamber_walls(cubic)
        # the spectral chain on a cheap Float64 harmonic well
        harm = SchrodingerProblem([0.0, 0.0, 1.0])
        spectral_cycles(harm, 1.0)
        wkb_eigenvalue(harm, 0, 0.3; order = 2)
        # the degenerate dual: a double turning point and its two refinements
        refinements(polygon_decomposition(stokes_graph(SchrodingerProblem([0.0, 0.0, 1.0]);
                                                       theta = 0.3)))
    end
end

end # module ExactWKB
