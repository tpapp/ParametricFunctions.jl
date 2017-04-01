immutable TrivialModel{T}
    α::T
    β::T
end

@testset "Collocation" begin
    fam = Chebyshev(10)
    θ = fit(fam, x->x^2)
    model = TrivialModel(2.0, 3.0)
    residual_function(model::TrivialModel, f, x) = f(x-model.α)-model.β
    cres = CollocationResidual(model, fam, residual_function)
    r = cres(θ)
    @test r ≈ (points(fam) - model.α).^2-model.β
end

immutable TrivialModel2{T}
    α::T
end

@testset "Collocation solver" begin
    (m::TrivialModel2)(x) = m.α*x
    model0 = TrivialModel2(2.0)
    residualf(model::TrivialModel2, f, x) = f(x)-model(x)
    res = CollocationResidual(model0, Chebyshev(10), residualf)
    f_sol, o = solve_collocation(res, zero)
    @test norm([f_sol(x)-model0(x) for x in linspace(domain(f_sol), 100)]) ≤ 1e-14
end
