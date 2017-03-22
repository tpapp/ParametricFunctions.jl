module ParametricFunctions

using ArgCheck
using ContinuousTransformations
import ContinuousTransformations: domain
using Lazy

export points, degf, basis, basis!, evaluate, basis_matrix, fit, fit!

abstract FunctionFamily

abstract ParametricFunction

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
basis{T}(p::ParametricFunction, x::T) = basis!(p, x, Vector{T}(degf(p)))

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

"""
Fit a parametric function to a function family, returning the
coefficients. See also the method for `\`.
"""
function fit{T}(p::ParametricFunction, y::AbstractVector{T})
    fit!(p, y, Vector{T}(degf(p)))
end

include("Chebyshev.jl")

end # module
