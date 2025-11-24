module GregoryLoredo


using CairoMakie
#using Format
using ProgressMeter
using Trapz


export FrequencyRange
export LightCurveShape
export MarginalizedPeriodogram
export OddRatios
export OddRatioSummary
export Periodogram
export PeriodogramSummary
export PlotLightCurve
export PlotOddRatios
export PlotPeriodogram



"""
  BinMult(nbins::Vector{Int})
  
#Arguments

- `nbins` is a vector with the number of events for each bin.

Compute the basic factor for the multiplicity of the binned events (Eq. 1.4 in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).
"""
function BinMult(nbins::Vector{Int})
  return sum(factorial.(BigInt.(nbins)))
end



"""
  BinningFactor(ω::Float64,Φ::Float64,m::Int,tlist::Vector{Float64})::BigFloat

#Arguments

- `ω` is the angular frequency.
- `Φ` is the phase.
- `m` is the number of bins.
- `tlist` is the vector of arrival times.

Compute the binning factor (Eq. B4 in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).
"""
function BinningFactor(ω::Float64,Φ::Float64,m::Int,tlist::Vector{Float64})::BigFloat
    njbins = ComputeBinValues(m,ω,Φ,tlist)
    return prod(WeightingFactor(njbins,m,length(tlist)).^-njbins)
end



"""
  ComputeBinValues(m::Int,ω::Float64,Φ::Float64,tlist::Vector{Float64})
  
#Arguments

- `m` is the number of bins.
- `ω` is the angular frequency.
- `Φ` is the phase.
- `tlist` is the vector of arrival times.

Compute the number of events falling in each bin of a model with `m` bins with given frequency and phase.
"""
function ComputeBinValues(m::Int,ω::Float64,Φ::Float64,tlist::Vector{Float64})
  ff = Int.(floor.(m*mod.(ω*tlist.+Φ,2π)/2π))
  return [count(==(i),ff) for i in 0:m-1]
end




"""
  FrequencyRange(Tlist::Vector{Float64};tstart=minimum(Tlist),tend=maximum(Tlist),logspacing=false))
  
#Arguments

- `Tlist` is the vector of arrival times.
- `tstart` is the minimum possible arrival time (default to the minimum of Tlist)
- `tend` is the minimum possible arrival time (default to the maximum of Tlist)
- `logspacing` generates a log-spaced frquency range.

Compute the frequency range to be analysed accordiong to the advises in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).

"""
function FrequencyRange(Tlist::Vector{Float64};tstart=minimum(Tlist),tend=maximum(Tlist),logspacing=false)
  N = length(Tlist)
  Δt = tend-tstart
  #
  ωhi = 2π*N/Δt
  ωlo = 20π/Δt
  dω = π/Δt   # step size
  if logspacing
      ωr = collect(logrange(ωlo,ωhi,length=N÷2))
  else
      ωr = collect(range(start=ωlo,stop=ωhi,step=dω))
  end
  return ωr
end





"""
  FullBinMult(nbins::Vector{Int})
  
#Arguments

- `nbins` is a vector with the number of events for each bin.

Compute the full multiplicity of the binned events (Eq. 1.4 in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).
"""
function FullBinMult(nbins::Vector{Int})
  N = sum(nbins)
  #
  return factorial(BigInt(N)) /  prod(factorial.(BigInt.(nbins)))
end



"""
  FullBinMultLog(nbins::Vector{Int})
  
#Arguments

- `nbins` is a vector with the number of events for each bin.

Compute the full multiplicity of the binned events (Eq. 1.4 in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).

The computation is here carried out by the logarithm.
"""
function FullBinMultLog(nbins::Vector{Int})
  N = sum(nbins)
  #
  return exp( log(factorial(BigInt(N))) - sum(log.(factorial.(BigInt.(nbins)))))
end




"""
  jt(m::Int,ω::Float64,Φ::Float64,t::Float64)
  
#Arguments

- `m` is the number of bins.
- `ω` is the angular frequency.
- `Φ` is the phase.
- `t` is the time.

Compute the bin number at time `t`for a model with `m` bins with given frequency and phase.
"""
function jt(m::Int,ω::Float64,Φ::Float64,t::Float64)
  ff = Int.(1 + floor.(m*mod.(ω*t.+Φ,2π)/2π))
  return ff
end




"""
  lcshape

# Fields
- μlc: Vector of rate mean value for each bin
- σlc: Vector of standard deviation for the rate for each bin
"""
struct lcshape
  μlc::Vector{Float64}
  σlc::Vector{Float64}
end



"""
  LightCurveShape(ωr::Vector{Float64},m::Int,tlist::Vector{Float64};unc=false,progress=false,phasesteps=10)::lcshape
  
#Arguments

- `ωr` is the vector of the analyzed frequencies.
- `m` is the number of bins.
- `tlist` is the vector of arrival times.
- `unc` computes light-curve uncertainty if selected.
- `progress` shows progress bars if selected.

Compute the light curve shape given the number `m` of bins and the input arrival time `tlist` (Eq. 7.10 in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).
"""
function LightCurveShape(ωr::Vector{Float64},m::Int,tlist::Vector{Float64};unc=false,progress=false,phasesteps=10)::lcshape
    intg = Vector{BigFloat}(undef,length(ωr))
    inintg = Vector{BigFloat}(undef,length(ωr))
    outintg = Vector{BigFloat}(undef,m)
    soutintg = Vector{BigFloat}(undef,m)
    #
    Φr = collect(range(start=0.,stop=2π/m,length=phasesteps))
    #
    p1 = Progress(length(ωr); enabled=progress)
    for ω in ωr
        intgs = [1/FullBinMult(ComputeBinValues(m,ω,Φ,tlist)) for Φ in Φr]
        integral = trapz(Φr, intgs)
        om = m * integral
        inintg[findfirst(x->x==ω, ωr)] = om/ω
        next!(p1)
    end
    integral = trapz(ωr, inintg)
    C = 1/integral  # Eq. 6.3 in Gregory & Lorendo (1992)
    #
    p2 = Progress(m; enabled=progress)
    for b in 1:m
        for ω in ωr
            intgs = [(1/ω) * (1/FullBinMult(ComputeBinValues(m,ω,Φ,tlist))) * LightCurveSimpleShape(ω,Φ,m,tlist).μlc[b] for Φ in Φr]
            integral = trapz(Φr, intgs)
            om = m * integral
            inintg[findfirst(x->x==ω, ωr)] = C*om
        end
        integral = trapz(ωr, inintg)
        outintg[b] = integral
        next!(p2)
    end
    #
    if unc == true
        p3 = Progress(m; enabled=progress)
        for b in 1:m
            for ω in ωr
                intgs = [(1/ω) * (1/FullBinMult(ComputeBinValues(m,ω,Φ,tlist))) * LightCurveSimpleShape(ω,Φ,m,tlist).μlc[b].^2 for Φ in Φr]
                integral = trapz(Φr, intgs)
                om = m * integral
                inintg[findfirst(x->x==ω, ωr)] = C*om
            end
            integral = trapz(ωr, inintg)
            soutintg[b] = integral
            next!(p3)
        end
    end
    #
    if unc == true
        return lcshape(outintg,sqrt.(soutintg .- outintg.^2))
    else
        return lcshape(outintg,zeros(length(outintg)))
    end
end





"""
  LightCurveSimpleShape(ω::Float64,Φ::Float64,m::Int,tlist::Vector{Float64};unc=false)::lcshape
  
#Arguments

- `ω` is the angular frequency.
- `Φ` is the phase.
- `m` is the number of bins.
- `tlist` is the vector of arrival times.
- `unc` computes light-curve uncertainty if selected.

Compute the light curve shape given the angular frequency `ω`, the phase `Φ`, `m` bins given the input arrival time `tlist` (Eq. 4.3 in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).
"""
function LightCurveSimpleShape(ω::Float64,Φ::Float64,m::Int,tlist::Vector{Float64};unc=false)::lcshape 
    μlc = ( ComputeBinValues(m,ω,Φ,tlist) .+ 1) ./ (length(tlist) + m)
    if unc == true
        σlc = sqrt.(μlc .* (1 .- μlc) ./ (length(tlist)+m+1))
    end
    #
    if unc == true
        return lcshape(μlc,σlc)
    else
        return lcshape(μlc,zeros(length(μlc)))
    end
end





"""
  logStirlingApprox(n::Int)

#Arguments

- `n` is the integer we want to use to compute the factorial.

"""
function logStirlingApprox(n::Int)
  if n == 0
    return 1
  else
    return n * log(n/ℯ) + 0.5*log(2π * n)
  end
end





"""
  MarginalizedPeriodogram(Tlist::Vector{Float64},ωr::Vector{Float64};progress=false,phasesteps=10)
  
#Arguments

- `Tlist` is the vector of arrival times.
- `ωr` is the vector of the analyzed frequencies.
- `mmax` is the maximum number of bins.
- `progress` shows progress bars if selected.


Compute the periodogram following Section 6 in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).

"""
function MarginalizedPeriodogram(Tlist::Vector{Float64},ωr::Vector{Float64},mmax::Int;progress=false,phasesteps=10)
    #
    intg = Vector{BigFloat}(undef,length(ωr))
    ppm = Array{Vector{BigFloat}}(undef,mmax-1,length(ωr))
    pppm = Vector{BigFloat}(undef,mmax-1)
    # Computing normalization factors... (eq. 6.3)
    Cm = []
    p1 = Progress(mmax-1; enabled=progress)
    for m in 2:mmax    
        Φr = collect(range(start=0.,stop=2π/m,length=phasesteps))
        for ω in ωr
            intgs = [1/FullBinMult(ComputeBinValues(m,ω,Φ,Tlist)) for Φ in Φr]
            integral = trapz(Φr, intgs)
            om = m * integral
            intg[findfirst(x->x==ω, ωr)] = om/ω
        end    
        integral = trapz(ωr, intg)
        C = 1/integral
        push!(Cm,C)
        next!(p1)
    end
    #
    # Computing p(ω)... (eq. 6.2)
    p2 = Progress(mmax-1; enabled=progress)
    for m in 2:mmax
        pp = Vector{BigFloat}(undef,m)  
        Φr = collect(range(start=0.,stop=2π/m,length=phasesteps))
        for ω in ωr
            intgs = [1/FullBinMult(ComputeBinValues(m,ω,Φ,Tlist)) for Φ in Φr]
            integral = trapz(Φr, intgs)
            intg[findfirst(x->x==ω, ωr)] = m * integral * Cm[m-1]
        end
        ppm[m-1] = intg
        next!(p2)
    end
    #
    N = length(Tlist)
    ωhi = maximum(ωr)
    ωlo = minimum(ωr)
    #
    omq = Vector{BigFloat}(undef,mmax-1)
    #
    factor1 = 1/(2π*(mmax-1))
    factor2 = 1/log(ωhi/ωlo)
    factor8 = factorial(BigInt(N))
    #
    # Computing Odds ratios... (Eq. 5.28)
    p3 = Progress(mmax-1; enabled=progress)
    for m in 2:mmax
        factor4 = factorial(BigInt(N+m-1))
        factor7 = exp(N*log(BigInt(m)))
        factor = factor1*factor2*factor8*factor7/factor4
        #
        Φr = collect(range(start=0.,stop=2π/m,length=phasesteps))
        M = [factor/FullBinMult(ComputeBinValues(m,ω,Φ,Tlist)) for Φ in Φr, ω in ωr]
        integral = trapz((Φr,ωr),M)
        om = m * integral
        omq[m-1] = om
        next!(p3)
    end
    # 
    # Computing model probabilities... (Eq. 2.13)
    sm = sum(omq) + 1.
    for m in 2:mmax
        pppm[m-1] = omq[m-1] / sm
    end
    #
    finper = zeros(BigFloat,length(ωr))
    for m in 2:mmax
        finper = finper .+ pppm[m-1] * ppm[m-1]
    end
    #
    return finper
end


"""
  orsum

# Fields
- prob: Vector of probabilities for each possible bin number `m`
- maxm: Bin number showing the maximum probability
"""
struct orsum
  prob::Float64
  maxm::Int
end



"""
  OddRatioSummary(or::Vector{Real})::orsum
  
#Arguments

- `od` is the vector of computed odd ratios for each considered light-curve bin number.

Compute the probability for a non constant model [Eq. 5.2. in Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)) and the bin number yileding the highest odd-ratio.

"""
function OddRatioSummary(or::Vector{BigFloat})::orsum
  return orsum(sum(or)/(1+sum(or)), argmax(or)+1)
end






"""
  OddRatios(Tlist::Vector{Float64},ωr::Vector{Float64},mmax::Int;progress=false,phasesteps=10)
  
#Arguments

- `Tlist` is the vector of arrival times.
- `ωr` is the vector of the analyzed frequencies.
- `mmax` is the maximum number of bins.
- `progress` shows progress bars if selected.


Compute the odd ratio following Appendix A in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).

"""
function OddRatios(Tlist::Vector{Float64},ωr::Vector{Float64},mmax::Int;progress=false,phasesteps=10)
  N = length(Tlist)
  ωhi = maximum(ωr)
  ωlo = minimum(ωr)
  #
  omq = Vector{BigFloat}(undef,mmax-1)
  intg = Vector{BigFloat}(undef,length(ωr))
  #
  factor1 = 1/(2π*(mmax-1))
  factor2 = 1/log(ωhi/ωlo)
  factor8 = factorial(BigInt(N))
  #
  p = Progress(mmax; enabled=progress)
  for m in 2:mmax
    #factor3 = factorial(BigInt(m-1))
    factor4 = factorial(BigInt(N+m-1))
    factor7 = exp(N*log(BigInt(m)))
    factor = factor1*factor8*factor7/factor4
    #
    Φr = collect(range(start=0.,stop=2π/m,length=phasesteps))
    for ω in ωr
        intgs = [factor/FullBinMult(ComputeBinValues(m,ω,Φ,Tlist)) for Φ in Φr]
        integral = trapz(Φr, intgs)
        om = m * integral
        intg[findfirst(x->x==ω, ωr)] = om * factor2 / ω
    end
    integral = trapz(ωr, intg)
    omq[m-1] = integral
    next!(p)
  end
  return omq
end



"""
  Periodogram(Tlist::Vector{Float64},ωr::Vector{Float64},m::Int,phasesteps=10)
  
#Arguments

- `Tlist` is the vector of arrival times.
- `ωr` is the vector of the analyzed frequencies.
- `m` is the number of bins.

Compute the periodogram following Section 6 in [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract)).

"""
function Periodogram(Tlist::Vector{Float64},ωr::Vector{Float64},m::Int,phasesteps=10)
    inintg = Vector{BigFloat}(undef,length(ωr))
    outintg = Vector{BigFloat}(undef,length(ωr))
    Φr = collect(range(start=0.,stop=2π/m,length=phasesteps))
    #
    for ω in ωr
        intgs = [1/FullBinMult(ComputeBinValues(m,ω,Φ,Tlist)) for Φ in Φr]
        integral = trapz(Φr, intgs)
        om = m * integral
        inintg[findfirst(x->x==ω, ωr)] = om/ω
        outintg[findfirst(x->x==ω, ωr)] = om
    end
    integral = trapz(ωr, inintg)
    C = 1/integral  # Eq. 6.3 in Gregory & Lorendo (1992)
    #
    for ω in ωr
        om = (C/ω) * m * outintg[findfirst(x->x==ω, ωr)]
        inintg[findfirst(x->x==ω, ωr)] = om
    end
    #
    return inintg
end


"""
  persum

# Fields
- maxpow: Maximum power in the periodogram
- maxfreq: Frequency at the maximum power
- maxper: Period at the maximum power
"""
struct persum
  maxpow::Float64
  maxfreq::Float64
  maxper::Float64
end



"""
  PeriodogramSummary(per::Vector{BigFloat},ωr::Vector{Float64})::persum
  
#Arguments

- `per` is the vector of powers for each considered frequency.
- `ωr` is the vector of the analyzed frequencies.


Plot the Odd Ratio as a function of the number of bins.

"""
function PeriodogramSummary(per::Vector{BigFloat},ωr::Vector{Float64})::persum
  return persum(per[argmax(per)],ωr[argmax(per)],2π/ωr[argmax(per)])
end




"""
  PlotLightCurve(lpl::lcshape;xlabel="",ylabel="",sigma=1.)
  
#Arguments

- `lpl` is a light-curve shape structure.
- `xlabel` is the (optional) label for the x-axis.
- `ylabel` is the (optional) label for the y-axis.
- `sigma` is the number od standard deviation to be plotted.



Plot the light-curve in the given bins of the analysis.



"""
function PlotLightCurve(lpl::lcshape;xlabel="",ylabel="",sigma=1.)
    f = Figure()

    μ = lpl.μlc
    σ = lpl.σlc
    
    ax = Axis(f[1,1],
    spinewidth=2,
    ylabel=ifelse(ylabel=="","Light Curve",ylabel),
    xlabel=ifelse(xlabel=="","Phase",xlabel),
    xlabelsize=20,
    ylabelsize=20,
    xgridvisible = false,
    ygridvisible = false,
    xticklabelsize = 20,
    yticklabelsize = 20,
    )

    p = stairs!((1:length(μ ))/length(μ ),μ )
    if count(x->x>0.,σ) > 0
      stairs!((1:length(μ ))/length(μ ),μ .- sigma*σ,linestyle=:dash, color=p.attributes.color)
      stairs!((1:length(μ ))/length(μ ),μ .+ sigma*σ,linestyle=:dash, color=p.attributes.color)
    end
    
    return f
end




"""
  PlotOddRatios(or::Vector{Real};xlabel="",ylabel="";xlabel="",ylabel="")
  
#Arguments

- `od` is the vector of computed odd ratios for each considered light-curve bin number.
- `xlabel` is the (optional) label for the x-axis.
- `ylabel` is the (optional) label for the y-axis.


Plot the oddratio for each possible bin number.

"""
function PlotOddRatios(or::Vector{BigFloat};xlabel="",ylabel="")
    f = Figure()
    
    ax = Axis(f[1,1],
    spinewidth=2,
    ylabel=ifelse(ylabel=="","Log Odd Ratio",ylabel),
    xlabel=ifelse(xlabel=="","Bin number",xlabel),
    xlabelsize=20,
    ylabelsize=20,
    xgridvisible = false,
    ygridvisible = false,
    xticklabelsize = 20,
    yticklabelsize = 20,
    xticks = 2:length(or)+1,
    )

    barplot!(2:length(or)+1,log10.(or))

    return f
end





"""
  PlotPeriodogram(per::Vector{BigFloat},ωr::Vector{Float64};xlabel="",ylabel="",pmax=2π/minimum(ωr),pmin=2π/maximum(ωr))
  
#Arguments

- `per` is the vector of powers for each considered frequency.
- `ωr` is the vector of the analyzed frequencies.
- `xlabel` is the (optional) label for the x-axis.
- `ylabel` is the (optional) label for the y-axis.
- `pmin` is the minimum of the plotted periods (default minimum of the input data).
- `pmax` is the maximum of the plotted periods (default maximum of the input data).


Plot the periodogram as a function of the analysed frequencies.



"""
function PlotPeriodogram(per::Vector{BigFloat},ωr::Vector{Float64};xlabel="",ylabel="",pmax=2π/minimum(ωr),pmin=2π/maximum(ωr))
    f = Figure()
    
    ax = Axis(f[1,1],
    spinewidth=2,
    ylabel=ifelse(ylabel=="","Log Periodogram",ylabel),
    xlabel=ifelse(xlabel=="","Period",xlabel),
    xlabelsize=20,
    ylabelsize=20,
    xgridvisible = false,
    ygridvisible = false,
    xticklabelsize = 20,
    yticklabelsize = 20,
    )
    
    ωrs = ωr[(ωr .> (2π / pmax)) .&& (ωr .< (2π / pmin))]
    pers = per[(ωr .> (2π / pmax)) .&& (ωr .< (2π / pmin))]

    lines!(2π ./ ωrs,log10.(pers))
    
    #xlims!(pmin,pmax)

    return f
end


"""
  WeightingFactor(njbin::Vector{Int},m::Int,totcnt::Int)::Vector{BigFloat}

#Arguments

- `njbin` is the vector with the number of events in each model bin.
- `m` is the number of bins for the adopted model.
- `totcnt` is the total number of events in the analyzed dataset.


Compute the weighting factors for each bin of the model.

"""
function WeightingFactor(njbin::Vector{Int},m::Int,totcnt::Int)::Vector{BigFloat}
    return ((njbin * m) / totcnt).^-njbin
end



end
