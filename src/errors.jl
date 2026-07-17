# Typed error hierarchy. Every invalid operation throws one of these; messages start
# with the error-type name so failures are self-identifying in test output. Generic
# bad-keyword cases reuse `Resurgence.InvalidArgument`.

"""
    ExactWKBError

Abstract supertype of all typed errors thrown by ExactWKB.jl.
"""
abstract type ExactWKBError <: Exception end

"""
    InvalidPotential(msg)

Thrown when a [`SchrodingerProblem`](@ref) is built from an inadmissible potential —
an empty coefficient vector, or a `Q = V − E` that is identically constant (no
turning points, no WKB problem).
"""
struct InvalidPotential <: ExactWKBError
    msg::String
end

Base.showerror(io::IO, e::InvalidPotential) = print(io, "InvalidPotential: ", e.msg)

"""
    UnsupportedTurningPoint(z, order)

Thrown when an operation that only handles simple turning points meets one of higher
`order` at `z` — e.g. the WKB recursion evaluation or the Stokes tracer at a double
turning point. The point is still classified correctly by [`turning_points`](@ref);
only the downstream local model is missing (deferred to M4).
"""
struct UnsupportedTurningPoint <: ExactWKBError
    z::Number
    order::Int
end

function Base.showerror(io::IO, e::UnsupportedTurningPoint)
    print(io, "UnsupportedTurningPoint: turning point at z ≈ $(e.z) has order ",
          "$(e.order); only simple (order-1) turning points are supported in M3")
end

"""
    ContourError(msg)

Thrown by period integration when a contour is ill-posed for the branch bookkeeping:
a vertex sitting on a turning point (for a closed contour), or a step across which the
√Q branch cannot be tracked unambiguously.
"""
struct ContourError <: ExactWKBError
    msg::String
end

Base.showerror(io::IO, e::ContourError) = print(io, "ContourError: ", e.msg)

"""
    TracingFailed(msg)

Thrown by the Stokes-graph tracer when a line neither hits a turning point nor
escapes to infinity within the mass budget (`max_mass`). Pass `allow_incomplete =
true` to [`stokes_graph`](@ref) to retain the partial line instead.
"""
struct TracingFailed <: ExactWKBError
    msg::String
end

Base.showerror(io::IO, e::TracingFailed) = print(io, "TracingFailed: ", e.msg)

"""
    NonGenericGraph(msg)

Thrown by [`ideal_triangulation`](@ref) when the traced Stokes graph is not generic —
it contains a finite Stokes line (`θ` sits on a wall: perturb it away from the
[`saddle_candidates`](@ref) phases), an `:incomplete` line, or fails a combinatorial
invariant of a generic degree-`d` graph (`d + 2` marked points, `d` triangles,
`2d + 1` Stokes regions).
"""
struct NonGenericGraph <: ExactWKBError
    msg::String
end

Base.showerror(io::IO, e::NonGenericGraph) = print(io, "NonGenericGraph: ", e.msg)

"""
    CoalescentTurningPoints(energy, z, separation)

Thrown by [`spectral_cycles`](@ref) when the energy `E` sits at (or numerically too
close to) a critical value of the potential, so two turning points coalesce near `z`
(separation `separation`; a multiple root reports `separation = 0`). The WKB cycle
data is singular there — perturb `E` away from the critical value. Degenerate-Weber
local models are deliberately out of the M5 scope.
"""
struct CoalescentTurningPoints <: ExactWKBError
    energy::Number
    z::Number
    separation::Real
end

function Base.showerror(io::IO, e::CoalescentTurningPoints)
    print(io, "CoalescentTurningPoints: at E = $(e.energy) two turning points ",
          "coalesce near z ≈ $(e.z) (separation ≈ $(e.separation)); E is too close ",
          "to a critical value of the potential — perturb the energy")
end

"""
    QuantizationError(msg)

Thrown by [`wkb_eigenvalue`](@ref) when the Newton iteration on the exact
quantization condition fails to converge (bad seed, coalescing turning points along
the Newton path, or an `n`/`parity` that selects no eigenvalue), and by
[`quantization_condition`](@ref) when the spectral-cycle layout at `E` does not match
the requested condition (no well cycle, or a multi-well problem where a single-well
condition was requested).
"""
struct QuantizationError <: ExactWKBError
    msg::String
end

Base.showerror(io::IO, e::QuantizationError) = print(io, "QuantizationError: ", e.msg)

"""
    ChamberError(msg)

Thrown when a chamber cannot carry the M4 cluster machinery: by
[`charge_basis`](@ref) for a *cyclic chamber* (a triangle of the triangulation with
all three edges internal, where the decay-oriented basis is not a seed basis — pick a
θ in a tree chamber instead), and by [`bps_spectrum`](@ref) when no valid chamber
ordering exists (two saddle phases coincide within tolerance — a marginal-stability
wall, perturb the potential or the energy — or the greedy θ-decreasing mutation
sequence fails to close into a maximal green sequence).
"""
struct ChamberError <: ExactWKBError
    msg::String
end

Base.showerror(io::IO, e::ChamberError) = print(io, "ChamberError: ", e.msg)
