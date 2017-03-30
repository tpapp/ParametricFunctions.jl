using RecipesBase
using Parameters

"""
Plot a univariate parametric function, marking the collocation points.

`N` is the number of points for the dense grid on which the function
is calculated. When and endpoint is infinite, `padding` is used to add
some space to the nearest collocation point, specified as a fraction
of the range covered by all collocation points. `label` labels the plot.
"""
@recipe function f(pf::ParametricFunction, label = "function", N = 100, padding = 0.1)
    @unpack family, θ = pf
    dom = domain(family)
    @assert isa(dom, AbstractInterval) "Plotting not implemented for domain $(dom)."
    left, right = extrema(dom)
    zs = points(family)
    if !isfinite(dom)
        pleft, pright = extrema(zs)
        pad = (pright - pleft)*padding
        if !isfinite(left)
            left = pleft - pad
        end
        if !isfinite(right)
            right = pright + pad
        end
    end
    xs = linspace(left, right, N)
    @series begin
        seriestype := :path
        label := label
        xs, pf.(xs)
    end
    @series begin
        primary := false
        seriestype := :scatter
        zs, pf.(zs)
    end
end
