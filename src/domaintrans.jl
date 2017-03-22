export DomainTrans

"""
Transform the domain of a `ParametricFunction`.

Evaluation at `x` with parameters `θ` is equivalent to
`evaluate(inner, θ, transformation(x))`.
"""
immutable DomainTrans{T, S} <: ParametricFamily
    transformation::T
    inner::S
end

function DomainTrans(dom::AbstractInterval, fam::ParametricFamily)
    DomainTrans(bridge(dom, domain(fam)), fam)
end

function DomainTrans(dom::AbstractInterval, transformation::UnivariateTransformation,
                     fam::ParametricFamily)
    DomainTrans(bridge(dom, transformation, domain(fam)), fam)
end

@forward DomainTrans.inner degf

domain(p::DomainTrans) = domain_in_image(p.transformation, domain(p.inner))

points(p::DomainTrans) = inv(p.transformation).(points(p.inner))

function basis!{T}(p::DomainTrans, x::T, b::AbstractVector{T})
    basis!(p.inner, p.transformation(x), b)
end

evaluate(p::DomainTrans, θ, x) = evaluate(p.inner, θ, p.transformation(x))

fit!(p::DomainTrans, ys, θ) = fit!(p.inner, ys, θ)
