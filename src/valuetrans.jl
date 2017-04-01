export ValueTrans

"""
Values obtained from the `inner` family are transformed by
`trans_function(x)`, which should be invertible when used with `fit`.

Basis is not defined, since the transformation is generally not
linear.
"""
immutable ValueTrans{T, S} <: ParametricFamily
    trans_function::T
    inner::S
end

@forward ValueTrans.inner degf, domain, points

function basis!(::ValueTrans, ::Any, ::Any)
    throw(ArgumentError("Value transformations are generally not linear, basis not defined."))
end

evaluate(p::ValueTrans, θ, x) = p.trans_function(x)(evaluate(p.inner, θ, x))

function fit!(p::ValueTrans, ys::AbstractVector, θ)
    zs = [inv(p.trans_function(x))(y) for (x,y) in zip(points(p.inner), ys)]
    fit!(p.inner, zs, θ)
end
