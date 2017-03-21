module ParametricFunctions

using ContinuousTransformations
import ContinuousTransformations: domain

using ArgCheck

export points, degf, basis, basis!, evaluate, basis_matrix, fit, fit!

abstract FunctionFamily

abstract ParametricFunction

"Domain of the parametric function or a function family."
function domain end

"Degrees of freedom (number of coefficients)."
function degf end

"Recommended points for collocation and function fitting."
function points end

function basis end

function basis! end

function evaluate end

function basis_matrix{T}(family, xs::AbstractVector{T})
    B = Array{T}(degf(family), length(xs))
    for (i, x) in enumerate(xs)
        basis!(family, x, @view B[i, :])
    end
    B
end

function fit! end

function fit end

include("Chebyshev.jl")

end # module
