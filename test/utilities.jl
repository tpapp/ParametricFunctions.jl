"Robust approximate comparison operator for two reals."
function ≅(a::Real, b::Real; ϵscale=1000, ϵpow=0.5)
    abs(a-b) < ϵscale*(max(eps(a),eps(b))^ϵpow)*(1+abs(a)+abs(b))
end

"""
Draw `n` random `θ` coefficients. Second parameter adds additional
dispersion on the log scale.
""" 
rand_θ(n, log_σ=4) = normalize!(randn(n) .* [exp(log_σ*randn()) for _ in n], 2)

"""
Accuracy of various f(θ, x), comparing the original and BigFloat
values at various `xs`. Report relative accuracy in digits, as quantiles q
"""
function accuracy(f, θ, xs, q=[0.0, 0.1, 0.5])
    y_bigf = f.([BigFloat.(θ)], BigFloat.(xs))
    y = f.([θ], xs)         # FIXME remove [] in v0.6
    acc(a,b) = -Float64(log10(abs(a-b)/(1+abs(Float64(a)))))
    quantile(acc.(y_bigf, y), q)
end

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
