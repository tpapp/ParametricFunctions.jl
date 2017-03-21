export Chebyshev

immutable Chebyshev
    n::Int
end

domain(::Chebyshev) = -1..1

points(c::Chebyshev, T = Float64) = cospi.(((c.n:-1:1)-T(0.5))/c.n)

degf(c::Chebyshev) = c.n

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

# note: after benchmarking, it was found that this is faster than Clenshaw.
# TODO: test for accuracy.
function evaluate{T}(c::Chebyshev, θ, x::T)
    xp = x
    xpp = one(T)
    s = zero(x)
    for i in 1:c.n
        if i == 1
            s += xpp*θ[i]
        elseif i == 2
            s += xp*θ[i]
        else
            xp, xpp = 2*x*xp - xpp, xp
            s += xp*θ[i]
        end
    end
    s
end

function fit!{T}(c::Chebyshev, y::AbstractVector{T}, θ::AbstractVector{T})
    @argcheck c.n == length(y) == length(θ)
    x = points(c, T)
    xpp = ones(T, c.n)
    xp = copy(x)
    for i in 1:c.n
        if i == 1
            θ[1] = sum(y) / c.n
        else
            if i > 2
                for j in 1:c.n
                    xp[j], xpp[j] = 2*x[j]*xp[j] - xpp[j], xp[j]
                end
            end
            θ[i] = dot(xp, y) *2 / c.n
        end
    end
    θ
end

fit{T}(c, y::AbstractVector{T}) = fit!(c, y, Vector{T}(c.n))

