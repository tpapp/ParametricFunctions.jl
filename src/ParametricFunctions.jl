module ParametricFunctions

using ContinuousTransformations
import ContinuousTransformations: domain

using ArgCheck

export points

abstract FunctionFamily

abstract ParametricFunction

"Domain of the parametric function or a function family."
function domain end

"Recommended points for collocation and function fitting."
function points end

include("Chebyshev.jl")

end # module
