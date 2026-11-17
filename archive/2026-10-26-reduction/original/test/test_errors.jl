@testset "errors" begin
    @test InvalidPotential <: ExactWKBError <: Exception
    @test UnsupportedTurningPoint <: ExactWKBError
    @test ContourError <: ExactWKBError
    @test TracingFailed <: ExactWKBError

    # showerror starts with the type name so failures self-identify
    for (e, tag) in ((InvalidPotential("x"), "InvalidPotential"),
                     (UnsupportedTurningPoint(0.0 + 0im, 2), "UnsupportedTurningPoint"),
                     (ContourError("y"), "ContourError"),
                     (TracingFailed("z"), "TracingFailed"))
        @test startswith(sprint(showerror, e), tag)
    end
    @test occursin("order 2", sprint(showerror, UnsupportedTurningPoint(1.0 + 0im, 2)))
end
