## Julia test script to run on Travis
## workaround until Pkg3 can REQUIRE master versions of packages

## for these packages, use master
Pkg.clone("https://github.com/tpapp/ContinuousTransformations.jl.git");

if VERSION ≥ v"0.6-"
    ## workaround for https://github.com/JuliaPlots/Plots.jl/issues/753
    for pkg in ["Plots", "PlotThemes", "PlotUtils"]
        Pkg.add(pkg)
        Pkg.checkout(pkg)
    end
end

## original test script
Pkg.clone(pwd()); Pkg.build("ParametricFunctions");
Pkg.test("ParametricFunctions"; coverage=true)
