"""
Calculate `∑ₖ cₖ fₖ, k ∈ 1:n` where `fₖ₊₁ = α fₖ + β fₖ₋₁` and the
first two elements are given. Uses the reverse Clenshaw algorithm,
stability is assumed.
"""
@inline function clenshaw{T}(c::AbstractVector{T}, α::T, β::T, f1::T, f2::T)
    bnn = bn = zero(T)
    for k in length(c):-1:2
        bn, bnn = c[k] + α*bn + β*bnn, bn
    end
    f1*(c[1]+β*bnn) + f2*bn
end

function evaluate(c::Chebyshev, θ, x, ::Type{Val{:clenshaw}})
    @argcheck c.n == length(θ)
    clenshaw(θ, 2*x, -one(x), one(x), x)
end
