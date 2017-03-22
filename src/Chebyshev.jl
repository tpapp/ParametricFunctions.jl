export Chebyshev

immutable Chebyshev <: ParametricFamily
    n::Int
end

domain(::Chebyshev) = -1..1

points(c::Chebyshev, T = Float64) = cospi.(((c.n:-1:1)-T(0.5))/c.n)

degf(c::Chebyshev) = c.n

function basis!{T}(c::Chebyshev, x::T, b::AbstractVector{T})
    @argcheck c.n == length(b)
    xp = x
    xpp = one(T)
    for i in 1:c.n
        if i == 1
            b[i] = xpp
        elseif i == 2
            b[i] = xp
        else
            xp, xpp = 2*x*xp - xpp, xp
            b[i] = xp
        end
    end
    b
end

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

function fit!{T}(c::Chebyshev, ys::AbstractVector{T}, θ::AbstractVector{T})
    @argcheck c.n == length(ys) == length(θ)
    x = points(c, T)
    xpp = ones(T, c.n)
    xp = copy(x)
    for i in 1:c.n
        if i == 1
            θ[1] = sum(ys) / c.n
        else
            if i > 2
                for j in 1:c.n
                    xp[j], xpp[j] = 2*x[j]*xp[j] - xpp[j], xp[j]
                end
            end
            θ[i] = dot(xp, ys) *2 / c.n
        end
    end
    θ
end
