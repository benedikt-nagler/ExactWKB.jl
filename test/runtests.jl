using ExactWKB
using Test
using Aqua
import Resurgence
include("diagonalization.jl")
@testset "ExactWKB.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(ExactWKB; ambiguities = false, persistent_tasks = (tmax = 300,))
        Aqua.test_ambiguities(ExactWKB)
    end
    include("test_errors.jl")
    include("test_potentials.jl")
    include("test_rational_potentials.jl")
    include("test_turning_points.jl")
    include("test_wkb_recursion.jl")
    include("test_periods.jl")
    include("test_voros.jl")
    include("test_weber.jl")
    include("test_quantum_periods.jl")
    include("test_quantization.jl")
    include("test_stokes_graph.jl")
    include("test_saddles.jl")
    include("test_ddp.jl")
    include("test_show.jl")
    include("test_makie_ext.jl")
end
