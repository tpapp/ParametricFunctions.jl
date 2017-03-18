using ParametricFunctions
using Base.Test
using ContinuousTransformations

@testset "Chebyshev polynomials" begin
    c = Chebyshev(9)

    @test domain(c) == -1..1
    @test all(p in -1..1 for p in points(c))
    @test points(c)[5] === 0.0
    @test basis(Chebyshev(5), 0.5) == [1.0, 0.5, -0.5, -1.0, -0.5]
end

