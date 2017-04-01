@testset "Chebyshev polynomials" begin
    c = Chebyshev(5)
    @test points(c)[3] === 0.0
    test_univariate(c, 5, -1..1; xs = linspace(-1,1,10), f = x->x^2+9*x-7)
end

@testset "Domain transformations" begin
    test_univariate(DomainTrans(ℝ⁺, Chebyshev(10)), 10, ℝ⁺; xs = linspace(0,5,10), f = x -> 1/(1+x)^2)
    test_univariate(DomainTrans(ℝ, REALCIRCLE, Chebyshev(10)), 10, ℝ;
                    xs = linspace(0,5,10), f = x -> 1/(1+x^2))
end

@testset "Value transformations" begin
    fam = ValueTrans(x->Shift(-log(x+one(x))), Chebyshev(10))
    test_univariate(fam, 10, -1..1;
                    xs = linspace(-1,1,10),
                    test_basis = false,
                    f = x -> 1/(2+x)^2+log(x+1),
                    test_partial = false)
    @test_throws ArgumentError basis(fam, 0.0)
end

@testset "Parametric function semantics" begin
    fam = Chebyshev(10)
    θ = rand_θ(10)
    f = ParametricFunction(fam, θ)
    @test degf(f) == 10
    @test parameters(f) == θ
    @test points(f) == points(fam)
    @test domain(f) == domain(fam)
end
