## Julia test script to run on Travis
## workaround until Pkg3 can REQUIRE master versions of packages

## for these packages, use master
Pkg.clone("https://github.com/tpapp/ContinuousTransformations.jl.git");
## workaround for https://github.com/JuliaPlots/Plots.jl/issues/753
## and https://github.com/JuliaImages/Images.jl/issues/604
if VERSION ≥ v"0.6-"
    for pkg in ["Plots", "PlotThemes", "PlotUtils", "Images"]
        Pkg.add(pkg)
        Pkg.checkout(pkg)
    end
end

## original test script
Pkg.clone(pwd()); Pkg.build("ParametricFunctions");
Pkg.test("ParametricFunctions"; coverage=true)

## exit
quit()
