export Chebyshev

immutable Chebyshev <: ParametricFamily
    n::Int
end

domain(::Chebyshev) = -1..1

points(c::Chebyshev, T = Float64) = cospi.(((c.n:-1:1)-T(0.5))/c.n)

degf(c::Chebyshev) = c.n

"""
Recurrence relation for Chebyshe polynomials.
"""
@inline _chebyshev_recurrence(x, T₋₁, T₋₂) = 2*x*T₋₁ - T₋₂, T₋₁

macro _chebyshev_recurrence!(x, T₋₁, T₋₂)
    :($T₋₁, $T₋₂ = _chebyshev_recurrence($x, $T₋₁, $T₋₂))
end

function basis!{T}(c::Chebyshev, x::T, b::AbstractVector{T})
    @argcheck c.n == length(b)
    T₋₁ = x
    T₋₂ = one(T)
    for i in 1:c.n
        if i == 1
            b[i] = T₋₂
        elseif i == 2
            b[i] = T₋₁
        else
            @_chebyshev_recurrence! x T₋₁ T₋₂
            b[i] = T₋₁
        end
    end
    b
end

# note: after benchmarking, it was found that this is faster than Clenshaw.
# TODO: test for accuracy.
function evaluate{T}(c::Chebyshev, θ, x::T)
    T₋₁ = x
    T₋₂ = one(T)
    s = zero(x)
    for i in 1:c.n
        if i == 1
            s += T₋₂*θ[i]
        elseif i == 2
            s += T₋₁*θ[i]
        else
            @_chebyshev_recurrence! x T₋₁ T₋₂
            s += T₋₁*θ[i]
        end
    end
    s
end

function fit!{T}(c::Chebyshev, ys::AbstractVector{T}, θ::AbstractVector{T})
    @argcheck c.n == length(ys) == length(θ)
    x = points(c, T)
    T₋₂ = ones(T, c.n)
    T₋₁ = copy(x)
    for i in 1:c.n
        if i == 1
            θ[1] = sum(ys) / c.n
        else
            if i > 2
                for j in 1:c.n
                    @_chebyshev_recurrence! x[j] T₋₁[j] T₋₂[j]
                end
            end
            θ[i] = dot(T₋₁, ys) *2 / c.n
        end
    end
    θ
end
