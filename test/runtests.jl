using SciMLTesting
using SciMLWorkshop
using Test

run_tests(
    ;
    core = () -> @test(nameof(SciMLWorkshop) === :SciMLWorkshop),
    qa = (;
        env = joinpath(@__DIR__, "qa"),
        body = joinpath(@__DIR__, "qa", "qa.jl"),
    ),
)
