export Chebyshev, basis

immutable Chebyshev
    n::Int
end

domain(::Chebyshev) = -1..1

points(c::Chebyshev) = cospi.(((c.n:-1:1)-0.5)/c.n)

function basis!{T}(c::Chebyshev, x::T, p::AbstractVector{T})
    xp = x
    xpp = one(T)
    for i in 1:c.n
        if i == 1
            p[i] = xpp
        elseif i == 2
            p[i] = xp
        else
            xp, xpp = 2*x*xp + xpp, xp
            p[i] = xp
        end
    end
    p
end
            
basis{T}(c::Chebyshev, x::T) = basis!(c, x, Vector{T}(c.n))
