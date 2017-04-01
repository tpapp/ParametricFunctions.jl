export Chebyshev

immutable Chebyshev <: ParametricFamily
    n::Int
end

domain(::Chebyshev) = -1..1

points(c::Chebyshev, T = Float64) = cospi.(((c.n:-1:1)-T(0.5))/c.n)

degf(c::Chebyshev) = c.n

"Recurrence relation for Chebyshev polynomials."
@inline _chebyshev_recurrence(x, T₋₁, T₋₂) = 2*x*T₋₁ - T₋₂, T₋₁

"Initialize T₋₁, T₋₂ for the Chebyshev recurrence."
@inline _chebyshev_recurrence_init{T}(x::T) = x, one(T)

macro _chebyshev_recurrence!(x, T₋₁, T₋₂)
    T₋₁ = esc(T₋₁)
    T₋₂ = esc(T₋₂)
    :(($(T₋₁), $(T₋₂)) = _chebyshev_recurrence($(esc(x)), $(T₋₁), $(T₋₂)))
end

function basis!{T}(c::Chebyshev, x::T, b::AbstractVector{T})
    @argcheck c.n == length(b)
    T₋₁, T₋₂ = _chebyshev_recurrence_init(x)
    for i in 1:c.n
        if i == 1
            b[i] = T₋₂
        else
            if i > 2
                @_chebyshev_recurrence! x T₋₁ T₋₂
            end
            b[i] = T₋₁
        end
    end
    b
end

# note: after benchmarking, it was found that this is faster than Clenshaw,
# and at least as accurate.
function evaluate{T}(c::Chebyshev, θ, x::T)
    T₋₁, T₋₂ = _chebyshev_recurrence_init(x)
    value = zero(T)
    for i in 1:c.n
        if i == 1
            value += T₋₂*θ[i]
        else
            if i > 2
                @_chebyshev_recurrence! x T₋₁ T₋₂
            end
            value += T₋₁*θ[i]
        end
    end
    value
end

function evaluate{T}(c::Chebyshev, θ, vp::ValuePartial{T})
    ## if T(x) = 2xT₋₁(x) - T₋₂(x)
    ## then T′(x) = 2T₋₁(x) + 2xT′₋₁(x) - T′₋₂(x), with T′₂ = 1 and T′₁ = 0
    x = vp.x
    value = partial = zero(T)
    T₋₁ = x
    T₋₂ = T′₋₁ = one(T)
    T′₋₂ = zero(T)
    for i in 1:c.n
        if i == 1
            value += T₋₂*θ[i]
        else
            if i > 2
                @_chebyshev_recurrence! x T′₋₁ T′₋₂
                T′₋₁ += 2*T₋₁
                @_chebyshev_recurrence! x T₋₁ T₋₂
            end
            value += T₋₁*θ[i]
            partial += T′₋₁*θ[i]
        end
    end
    value, partial
end

evaluate{T}(c::Chebyshev, θ, x::Partial{T}) = evaluate(c, θ, ValuePartial(x.x))[2]

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
