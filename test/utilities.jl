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
