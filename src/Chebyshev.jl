export Chebyshev, basis

immutable Chebyshev
    n::Int
end

domain(::Chebyshev) = -1..1

points(c::Chebyshev) = cospi.(((c.n:-1:1)-0.5)/c.n)

function basis!{T}(c::Chebyshev, x::T, p::AbstractVector{T})
    @argcheck c.n == length(p)
    xp = x
    xpp = one(T)
    for i in 1:c.n
        if i == 1
            p[i] = xpp
        elseif i == 2
            p[i] = xp
        else
            xp, xpp = 2*x*xp - xpp, xp
            p[i] = xp
        end
    end
    p
end

basis{T}(c::Chebyshev, x::T) = basis!(c, x, Vector{T}(c.n))

"""
Calculate ``∑ₖ cₖ fₖ, k=0,...,n`` where ``fₖ₊₁ = α fₖ + β fₖ₋₁`` and
the first two elements are given. Uses the reverse Clenshaw algorithm,
stability is assumed.
"""
function clenshaw{T}(c::AbstractVector{T}, α::T, β::T, f1::T, f2::T)
    bnn = bn = zero(T)
    for k in length(c):-1:2
        bn, bnn = c[k] + α*bn + β*bnn, bn
    end
    f1*(c[1]+β*bnn) + f2*bn
end

function evaluate(c::Chebyshev, θ, x)
    @argcheck c.n == length(θ)
    clenshaw(c, 2*x, -one(x), one(x), x)
end
