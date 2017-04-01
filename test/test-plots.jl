using VisualRegressionTests
using Plots; gr()

@testset "Plots" begin
    ref_image = Pkg.dir("ParametricFunctions", "test", "img00.png")
    function test_plot(fn)
        plt = plot(fitfun(Chebyshev(10), exp), label = false)
        png(plt, fn)
    end
    # test_plot(ref_image) # regenerates
    @test test_images(VisualTest(test_plot, ref_image), popup = false) |> success
end
