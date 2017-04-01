using ParametricFunctions
using Base.Test
using ContinuousTransformations
import ForwardDiff: derivative

include("utilities.jl")

include("test-accuracy.jl")
include("test-families.jl")
include("test-collocation.jl")
include("test-plots.jl")
