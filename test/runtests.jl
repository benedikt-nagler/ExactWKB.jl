using ExactWKB
using Test
using Aqua
import Resurgence

@testset "ExactWKB.jl" begin
    @testset "Aqua" begin
        # ClusterAlgebras is a declared dep for the M4 bridge but unused in M3.
        Aqua.test_all(ExactWKB; ambiguities = false,
                      stale_deps = (ignore = [:ClusterAlgebras],))
        Aqua.test_ambiguities(ExactWKB)
    end
    include("test_errors.jl")
    include("test_potentials.jl")
    include("test_turning_points.jl")
    include("test_wkb_recursion.jl")
    include("test_periods.jl")
    include("test_voros.jl")
    include("test_stokes_graph.jl")
    include("test_saddles.jl")
    include("test_show.jl")
    include("test_makie_ext.jl")
end
