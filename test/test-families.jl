function test_univariate(fam, expected_degf, expected_domain;
                         xs = nothing, test_basis = true, f = nothing,
                         test_partial = true, test_valuepartial = test_partial)
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
        θ0 = zeros(degf(fam))
        fit!(fam, f, θ0)
        θ1 = fit(fam, f)

        @test θ0 == θ1
        
        for z in zs
            @test evaluate(fam, θ1, z) ≅ f(z)
        end
        
        if test_basis
            θ2 = basis_matrix(fam, zs) \ ys
            @test θ1 ≈ θ2
        end
    end

    # when f given, test function fitting
    if f ≠ nothing
        zs = points(fam)
        pf = fitfun(fam, f.(zs))

        @test points(pf) == points(fam)
        @test degf(pf) == degf(fam)
        @test domain(pf) == domain(fam)
        @test family(pf) == fam
        
        for z in zs
            fz = f(z)
            f′z = derivative(f, z)
            @test pf(z) ≅ fz
            if test_partial
                @test pf(Partial(z)) ≅ f′z
            end
            if test_valuepartial
                v, p = pf(ValuePartial(z))
                @test v ≅ fz
                @test p ≅ f′z
            end
        end
    end
end

@testset "Chebyshev polynomials" begin
    c = Chebyshev(5)
    @test points(c)[3] === 0.0
    test_univariate(c, 5, -1..1; xs = linspace(-1,1,10), f = x->x^2+9*x-7)
end

@testset "Domain transformations" begin
    test_univariate(DomainTrans(ℝ⁺, Chebyshev(10)), 10, ℝ⁺; xs = linspace(0,5,10), f = x -> 1/(1+x)^2)
    test_univariate(DomainTrans(ℝ, REALCIRCLE, Chebyshev(10)), 10, ℝ; xs = linspace(0,5,10), f = x -> 1/(1+x)^2)
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
