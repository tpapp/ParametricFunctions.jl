module ParametricFunctions

using ArgCheck
using ContinuousTransformations
import ContinuousTransformations: domain
using Lazy

export                          # FIXME reorder nicely when interface stabilizes
    points, degf, basis, basis!, Partial, ValuePartial, evaluate,
    basis_matrix, fit, fit!, ParametricFamily, family, parameters, fitfun,
    ParametricFunction

abstract FunctionFamily

abstract ParametricFamily

"Domain of the parametric function or a function family."
function domain end

"Degrees of freedom (number of coefficients)."
function degf end

"Recommended points for collocation and function fitting."
function points end

"""
Basis functions evaluated at a given point, written into the third
argument which should be a vector.
"""
function basis! end

"""
Return the basis functions of `p` evaluated at `x`, as a vector.
"""
basis{T}(p::ParametricFamily, x::T) = basis!(p, x, Vector{T}(degf(p)))

"""
Partial derivative in the given coordinate.
"""
immutable Partial{T} x::T end

"""
Value and partial derivative in the given coordinate.
"""
immutable ValuePartial{T} x::T end

"""
Evaluate a parametric function with the given parameters.
"""
function evaluate end

"""
Return the basis matrix of a parametric function family evaluated at
`xs`. Note that the returned value may or may not be a dense matrix,
but will always be a conformable AbstractMatrix.
"""
function basis_matrix{T}(family, xs::AbstractVector{T})
    B = Array{T}(degf(family), length(xs))
    for (i, x) in enumerate(xs)
        basis!(family, x, @view B[i, :])
    end
    B
end

"""
Fit a parametric function to a function family, placing the
coefficients in the third argument.
"""
function fit! end

"Fit a function to a parametric family."
fit!{T <: ParametricFamily}(p::T, f::Function, θ) =  fit!(p, f.(points(p)), θ)

"""
Fit a parametric function to a function family, returning the
coefficients.
"""
fit(p::ParametricFamily, ys::AbstractVector) = fit!(p, ys, Vector{eltype(ys)}(degf(p)))

fit(p::ParametricFamily, f::Function) = fit(p, f.(points(p)))

immutable ParametricFunction{TF <: ParametricFamily, Tθ <: AbstractVector}
    family::TF
    θ::Tθ
    function ParametricFunction(fam, θ)
        @argcheck degf(fam) == length(θ) "Incompatible coefficients."
        new(fam, θ)
    end
end

ParametricFunction{TF, Tθ}(fam::TF, θ::Tθ) = ParametricFunction{TF, Tθ}(fam, θ)

"Parametric family of a parametric function."
family(f::ParametricFunction) = f.family

"Parameters of a parametric function."
parameters(f::ParametricFunction) = f.θ

@forward ParametricFunction.family domain, points, degf

"""
Fit a function from a parametric family at values `ys`, which were
evaluated at `points(family)`.
"""
fitfun(family::ParametricFamily, ys) = ParametricFunction(family, fit(family, ys))

(f::ParametricFunction)(x) = evaluate(f.family, f.θ, x)

include("Chebyshev.jl")
include("domaintrans.jl")
include("valuetrans.jl")
include("collocation.jl")
include("plots.jl")

end # module
