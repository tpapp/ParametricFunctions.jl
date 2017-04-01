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
