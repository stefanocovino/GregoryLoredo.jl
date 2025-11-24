using GregoryLoredo
using CairoMakie
using Test

@testset "GregoryLoredo.jl" begin
    # Write your tests here.
    @test GregoryLoredo.BinMult([1,4,6,10,4,2]) == 3629571
    #
    tlist1 = [0.3,1.2,3.4,5.,7.8,10.2,12.4,15.3]
    @test isapprox(GregoryLoredo.BinningFactor(2π/3.5,0.25,3,tlist1), 710.99, rtol=1e-4)
    #
    @test GregoryLoredo.ComputeBinValues(3,2π/3.5,0.25,tlist1) == [3,4,1]
    #
    tlist2 = [0.3,1.2,3.4,5.,7.8,10.2,12.4,15.3,18.5,20.3,23.2,27.3]
    @test FrequencyRange(tlist2) == [ 2.3271056693257726, 2.443460952792061, 2.55981623625835, 2.6761715197246385, 2.792526803190927]
    #
    @test GregoryLoredo.FullBinMult(GregoryLoredo.ComputeBinValues(3,2.44,0.25,tlist2)) == 27720.0
    #
    @test isapprox(GregoryLoredo.FullBinMultLog(GregoryLoredo.ComputeBinValues(3,2.44,0.25,tlist2)),27720.0, rtol=1e-8)
    #
    @test GregoryLoredo.jt(3,2.44,0.25,3.4) == 2
    #
    fr = FrequencyRange(tlist2)
    lcs = LightCurveShape(fr,3,tlist2)
    @test lcs.μlc == [0.4520096477412845, 0.3894071173763123, 0.15858323488240322]
    #
    lcss = GregoryLoredo.LightCurveSimpleShape(2.44,0.25,3,tlist2)
    @test lcss.μlc == [0.4, 0.26666666666666666, 0.3333333333333333]
    #
    @test GregoryLoredo.logStirlingApprox(10) == 15.096082009642153
    #
    @test isapprox(MarginalizedPeriodogram((tlist2),fr,4), [0.6138367,2.166482625,11.7351142,0.541620656,1.474411786574], rtol=1e-5)
    #
    or = OddRatios(tlist2,fr,4)
    @test isapprox(or, [0.5814103145,0.2704531661,0.0880807], rtol=1e-5)
    #
    pr = Periodogram(tlist2,fr,4)
    @test isapprox(pr, [1.4940811,5.022121,25.9666502,1.1463537898,2.9906035], rtol=1e-5)
    #
    @test typeof(PlotLightCurve(lcs)) == Figure
    #
    @test typeof(PlotOddRatios(or)) == Figure
    #
    @test typeof(PlotPeriodogram(pr,fr)) == Figure
    #
end
