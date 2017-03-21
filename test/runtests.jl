using ParametricFunctions
using Base.Test
using ContinuousTransformations

@testset "Chebyshev polynomials" begin
    c = Chebyshev(5)

    @test domain(c) == -1..1
    @test all(p in -1..1 for p in points(c))
    @test points(c)[3] === 0.0
    @test basis(c, 0.5) == [1.0, 0.5, -0.5, -1.0, -0.5]
    θ =  [0.083995, 0.565214, 0.467572, 0.318846, 0.0289252]
    @test evaluate(c, θ, 0.5) == dot(basis(c, 0.5), θ)

    p = points(c)
    B = basis_matrix(c, p)
    f = exp
    θ1 = B \ f.(p)

    @test evaluate(c, θ1, 0.0) ≈ f(0.0)
    
    θ2 = fit(c, f.(p))

    @test θ1 ≈ θ2

    @test evaluate(c, θ2, 0.0) ≈ f(0.0)
   
end

