import NLsolve

export CollocationResidual, solve_collocation

"""
Return an object that, when given parameters `θ`, construct a
parametric function using that an the family `fam`, and evaluates
`residual(model, f, x)` at the collocation points `x` of
`fam`. The two-argument version puts the residuals in the second
argument (should be a vector), while the first one returns them.

NOTE: Useful for passing to a solver such as `NLsolve.nlsolve`.
"""
immutable CollocationResidual{Tmodel, Tfam <: ParametricFamily, Txs, Tresidual} <: Function
    model::Tmodel
    fam::Tfam
    xs::Txs
    residual::Tresidual
    function CollocationResidual(model, fam, xs, residual)
        @argcheck degf(fam) == length(xs)
        new(model, fam, xs, residual)
    end
end

function CollocationResidual{Tmodel, Tfam, Txs, Tresidual}(model::Tmodel, fam::Tfam,
                                                           xs::Txs, residual::Tresidual)
    CollocationResidual{Tmodel, Tfam, Txs, Tresidual}(model, fam, xs, residual)
end

function CollocationResidual(model, fam, residual)
    CollocationResidual(model, fam, points(fam), residual)
end

@forward CollocationResidual.fam degf

function (cr::CollocationResidual)(θ, r)
    f = ParametricFunction(cr.fam, θ)
    r .= cr.residual.([cr.model], [f], cr.xs)
end

(cr::CollocationResidual)(θ) = cr(θ, Vector{eltype(cr.xs)}(degf(cr)))

function solve_collocation(cr::CollocationResidual, f₀; autodiff = true, options...)
    θ₀ = fit(cr.fam, f₀)
    o = NLsolve.nlsolve(cr, θ₀, autodiff = autodiff, options...)
    !NLsolve.converged(o) || warn("optimization did not converge")
    ParametricFunction(cr.fam, o.zero), o
end
