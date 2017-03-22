using ParametricFunctions
using Base.Test
using ContinuousTransformations

function test_univariate(fam, expected_degf, expected_domain;
                         xs = nothing, test_basis = true, f = nothing)
    @test degf(fam) == expected_degf
    @test domain(fam) == expected_domain
    @test all(p ∈ expected_domain for p in points(fam))

    # when xs given, test evaluation using basis
    if xs ≠ nothing
        θ = rand(expected_degf)
        for x in xs
            if test_basis
                @test evaluate(fam, θ, x) == dot(basis(fam, x), θ)
            else
                @test isa(evaluate(fam, θ, x), typeof(x))
            end
        end
    end

    # when f given, test approximation
    if f ≠ nothing
        zs = points(fam)
        ys = f.(zs)
        θ1 = fit(fam, f.(zs))
        
        for z in zs
            @test evaluate(fam, θ1, z) ≈ f(z)
        end
        
        if test_basis
            θ2 = basis_matrix(fam, zs) \ ys
            @test θ1 ≈ θ2
        end
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

@testset "Value transformations" begin
    fam = ValueTrans(x->Shift(-log(x+one(x))), Chebyshev(10))
    test_univariate(fam, 10, -1..1;
                    xs = linspace(-1,1,10),
                    test_basis = false,
                    f = x -> 1/(2+x)^2+log(x+1))
    @test_throws ArgumentError basis(fam, 0.0)
end
