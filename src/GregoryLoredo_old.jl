module GregoryLoredo


using Trapz


export compute_bin
export compute_GL



"""
  compute_bin(Tlist::Vector{Float64},m::Int,w::Float64,p::Float64=0.)

#Arguments

- `Tlist` is the vector of arrival times.
- `m` is the (optimal) bin size.
- `w` is the (peak) frequency
- `p` is the phase (default 0.).


"""
function compute_bin(Tlist::Vector{Float64},m::Int,w::Float64,p::Float64=0.)
  n = zeros(Int,m)
  j = floor.(m*mod.(w*Tlist.+p, 2*pi) / (2*pi))
  ji = Int.(j)
  for u in 1:m
    n[u] = length(ji[ji.==u-1])
  end
  return n
end




"""
  compute_factor(N::Int,m::Int,v::Float64)

#Arguments

- `N` is the number of input points.
- `m` is the bin number
- `v` is the maximum number of bins.


Compute ``m^N * (N+m-1)`` over ``N /(2 π v)`` which is used to scale the multiplicity function;
see Eq.5.25 of the original [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract) paper.

We return the log of the integral scale here to avoid numerical overflow.

"""
function compute_factor(N::Int,m::Int,v::Int)
  f1 = N*log(m)            # log of m^N
  f2 = sum(log.(1:N+m-1))  # log of (N+m-1)!
  f3 = sum(log.(1:m+1-1))  # log of m!
  f=f1+f3-f2-log(2*pi*v)
  return f
end



"""
  compute_GL(Tlist::Vector{Float64};w_range=nothing,m_max=12,ni=10)

#Arguments

- `Tlist` is the vector of arrival times.
- `w_range` is the vector of frequencies.
- `m_max` is the maximum number of bins (default=12).
- `ni` is the number of integration steps (default=10).


Compute the Gregory-Laredo algorithm on arrival times

This function computes the likelihood of a set of arrival times originating
from a periodic system rather than constant rate (poisson) background noise
based on [Gregory & Lorendo (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract).


Output:
O_period - Odds ratio for a periodic process vs. constant rate process
p_period - probability of a periodic process 0<=p_period<=1
m_opt    - optimal bin size 1<= m_opt <=m_max
S        - The probability spectrum
w        - The frequency range for S
w_peak   - the expected frequency of the process
w_conf   - 95% confidence interval of w_peak
"""
function compute_GL(Tlist::Vector{Float64};w_range=nothing,m_max=12,ni=10)
  N=length(Tlist)
  if N > 0
    fbin = precompute_binmult(N)
    v = m_max-1
    T= maximum(Tlist)-minimum(Tlist)
    if w_range === nothing               # use default frequencies
      w_hi = pi * N / T                 # max default frequency
      w_lo = minimum([20,N/10])*pi/T    # min default frequency
      dw = pi/T/10                      # step size
      w = collect(range(start=w_lo,stop=w_hi,step=dw))[1:end-1]
      if length(w) < 2
        error("Bad arrival time list.")
      end
    else                                # use user supplied frequency vector
      w = w_range
      w_hi = maximum(w_range)
      w_lo = minimum(w_range)
    end
    if w_lo == 0                        # cannot have w_lo =0
      error("Minimum frequency cannot be 0!\n")
    end
    #
    fa = zeros(m_max)
    for m in range(1,m_max)             # precompute factors for each m
      fa[m] = compute_factor(N,m,v)
    end
    #
    Om1w = compute_Om1w(Tlist,m_max,w,fa,fbin,ni)
    #println(Om1w)
    #
    pw = 1 ./ w ./ log(w_hi / w_lo)
    O1m = zeros(m_max)
    for i in 1:m_max                    # integrate odd ratios across frequencies
        O1m[i] =  trapz(w,pw.*Om1w[i,:])
    end
    #
    m_opt = argmax(O1m)                # find optimum bin number, i.e largest odds-ratio
    S = Om1w[m_opt,:] ./ w               # compute Spectral probability
    #m_opt = m_opt + 1                 # start bin index with 1
    C = trapz(w,S)                      # compute normalization
    S = S/C                             # normalized probability
    O_period = sum(O1m)                 # integrated odds ratio
    p_period = O_period / (1+O_period)  # likelihood of periodic event
    cdf = copy(S)
    for i in 1:length(S)
      cdf[i]= trapz(w[1:i-1],S[1:i-1])
    end
    wr = w[(cdf .> 0.025) .&& (cdf .< 0.975)]
    w_peak = w[argmax(S)]
    w_mean = trapz(w,S.*w)
    if length(wr)>0
      w_conf = [minimum(wr),maximum(wr)]
    else
      w_conf = [w_peak,w_peak]
    end
    return O_period,p_period,m_opt,S,w,w_peak,w_mean,w_conf
  else
    # throw an error
    error("No valid arrival time array provided!\n")
  end
end








"""
  compute_Om1(w::Float64,Tlist::Vector{Float64},m::Int,factor::Float64,fbin::Vector{Float64},ni::Int)
 
  
#Arguments

- `w` is the frequency.
- `Tlist` is the vector of arrival times.
- `m` is the maximum bin size.
- `factor` is the factor.
- `fbin` is the vector of bin factorial.
- `ni` is the number of integration steps.


Compute specific odds-ratios (eqn. 5.25).
"""
function compute_Om1(w::Float64,Tlist::Vector{Float64},m::Int,factor::Float64,fbin::Vector{Float64},ni::Int)
  p = 2*pi*collect(0:ni-1)./(ni*m)  # integration range, only integrate over a single bin, as values of the integral repeat over bins
  y = zeros(Float64,length(p))      # array to hold values of W_scaled over the integration range
  for i in 1:length(y)
    y[i] = GregoryLoredo.compute_W_scaled(Tlist,m,w,p[i],factor,fbin)
  end
  return trapz(p,y) * m             # return integrated W_Scaled
end





"""
  compute_Om1w(Tlist::Vector{Float64},m_max::Int,w::Vector{Float64},fa::Vector{Float64},fbin::Vector{Float64},ni::Int)

#Arguments

- `Tlist` is the vector of arrival times.
- `m_max` is the maximum number of bins.
- `w` is the vector of frequencies.
- `fa` is the vector of factors to scale the multiplicity function.
- `fbin` is the vector of factorial logarithms.
- `ni` is the number of integration steps.

Compute odds-ratios for bins and frequencies
"""
function compute_Om1w(Tlist::Vector{Float64},m_max::Int,w::Vector{Float64},fa::Vector{Float64},fbin::Vector{Float64},ni::Int)
    Om1w = zeros(Float64,(m_max,length(w))) # odds ratio matrix
    for m in 1:m_max
        for wi in 1:length(w)
            Om1w[m,wi] = compute_Om1(w[wi],Tlist,m,fa[m],fbin,ni)
        end
    end
    return Om1w
end
                    



	
"""
  compute_W_scaled(Tlist::Vector{Float64},m::Int,w::Float64,phase::Float64,factor::Float64,fbin::Vector{Float64})

#Arguments

- `Tlist` is the vector of arrival times.
- `m` is the maximum bin size.
- `w` is the frequency.
- `phase` is the phase for computation.
- `factor` is the factor.
- `fbin` is the vector of bin factorials.


Compute the scaled multiplicity 1/Wm(w,phase) (see eqn. 5.25)
  note that for large arrival time numbers the multiplicity becomes
  excessively large. Since W_m(phase,w) is never needed explicitly,
  we use the scaled version.
  input factor is the log of the actual factor
"""
function compute_W_scaled(Tlist::Vector{Float64},m::Int,w::Float64,phase::Float64,factor::Float64,fbin::Vector{Float64})
    n = compute_bin(Tlist,m,w,phase) # find bin histogram for a given bin number m, frequency w and phase
    f=0
    for i in 1:m
        f = f + fbin[n[i]+1]
    end
    return exp(f + factor)
end



"""
  precompute_binmul(N::int)
  
#Arguments

- `N` is the number of bins.

Precompute all potential bin factorials.
"""
function precompute_binmult(N::Int)
	#
	fbin = zeros(N+1)
	for i in range(start=3,stop=N+1) # n=0 -> define log(n)=0, n=1, log(n)=0, so no need to compute n=0,1
	  fbin[i] = fbin[i-1] + log(i-1)
	end
	return fbin
end






end