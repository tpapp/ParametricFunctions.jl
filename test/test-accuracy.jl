@testset "Chebyshev accuracy" begin
    for i in 1:100
        @test accuracy((θ,x)->evaluate(Chebyshev(length(θ)), θ, x), rand_θ(10),
                       linspace(-1,1,1001), [0.1])[1] ≥ 14
    end
    for i in 1:100
        @test accuracy((θ,x)->evaluate(Chebyshev(length(θ)), θ, Partial(x)), rand_θ(10),
                       linspace(-1,1,1001), [0.1])[1] ≥ 13
    end
end
