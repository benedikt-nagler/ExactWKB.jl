"""
    ExactWKB

Compute finite WKB expansions, numerical periods and selected Stokes graphs for
one-dimensional Schrödinger problems. Quantization is a bounded numerical approximation.
"""
module ExactWKB
using PrecompileTools: @setup_workload, @compile_workload
import Resurgence
import ClusterAlgebras

export ExactWKBError, InvalidPotential, UnsupportedTurningPoint, ContourError,
       TracingFailed, CoalescentTurningPoints, QuantizationError
export AbstractSchrodingerProblem, SchrodingerProblem, RationalProblem,
       q_coefficients, q_numerator, energy, degree, variable, with_energy,
       q_derivative_at, q_taylor_at, asymptotic_directions,
       poles, pole_orders, n_finite_poles, puncture_indices, n_punctures
export TurningPoint, turning_points, simple_turning_points, location, order, is_simple
export WKBExpansion, wkb_expansion, s_odd_terms, evaluate_s_odd, even_odd_residual
export WKBDerivative, wkb_derivative
export period_integral, wkb_period, encircling_contour, period_derivative
export VorosSymbol, voros_symbol, classical_period, quantum_series, full_series, voros_value
export SpectralCycle, spectral_cycles, quantum_period, kind, endpoints
export quantization_condition, spectral_determinant, wkb_eigenvalue
export weber_problem, weber_voros_series, weber_voros_sum, weber_log_gamma, weber_connection
export StokesLine, StokesGraph, stokes_graph, is_finite_line,
       finite_lines, edges, n_infinite_lines, topology_signature, is_degenerate,
       turning_point_orders, mass, endpoint
export Saddle, saddle_candidates, saddles, is_saddle, stokes_graph_family, central_charge
export ddp_transform, verify_ddp, verify_ddp_mutation, ddp_seed, intersection_pairing
export plot_stokes_graph

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
include("show.jl")

"""
    plot_stokes_graph(g; kwargs...)

Plot a Stokes graph with a loaded Makie backend.
"""
function plot_stokes_graph(args...; kwargs...)
    error("plot_stokes_graph requires a Makie backend; load CairoMakie or GLMakie")
end

@setup_workload begin
    @compile_workload begin
        p = SchrodingerProblem([-1.0, 0.0, 1.0])
        w = wkb_expansion(p; order = 2)
        voros_symbol(w, encircling_contour(turning_points(p)...))
        topology_signature(stokes_graph(SchrodingerProblem([0.0, 1.0])))
    end
end
end # module ExactWKB
