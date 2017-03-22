using ParametricFunctions
using Base.Test
using ContinuousTransformations

function test_univariate(fam, expected_degf, expected_domain;
                         xs = nothing, f = nothing, )
    @test degf(fam) == expected_degf
    @test domain(fam) == expected_domain
    @test all(p ∈ expected_domain for p in points(fam))

    # when xs given, test evaluation using basis
    if xs ≠ nothing
        θ = rand(expected_degf)
        for x in xs
            @test evaluate(fam, θ, x) == dot(basis(fam, x), θ)
        end
    end

    # when f given, test approximation
    if f ≠ nothing
        zs = points(fam)
        ys = f.(zs)
        B = basis_matrix(fam, zs)
        θ1 = B \ ys
        
        for z in zs
            @test evaluate(fam, θ1, z) ≈ f(z)
        end
        
        θ2 = fit(fam, f.(zs))
        
        @test θ1 ≈ θ2
    end
end

@testset "Chebyshev polynomials" begin
    c = Chebyshev(5)
    @test points(c)[3] === 0.0
    test_univariate(c, 5, -1..1; xs = linspace(-1,1,10), f = exp)
end

@testset "Domain transformations" begin
    fam = DomainTrans(ℝ⁺, Chebyshev(10))
    test_univariate(fam, 10, ℝ⁺; xs = linspace(0,5,10), f = x -> 1/(1+x)^2)
end
