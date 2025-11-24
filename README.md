# GregoryLoredo

[![Build Status](https://github.com/stefanocovino/GregoryLoredo.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stefanocovino/GregoryLoredo.jl/actions/workflows/CI.yml?query=branch%3Amain)


Phil Gregory and Tom Loredo published in 1992 a [seminal paper](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..146G/abstract) about a Bayesian method to detect periodic signal in astronomical datasets without assuming a particular shape for the periodic pattern.

The alghorithm has beed widely applied in astronomical researche and this is its implementation in the `Julia` language.

An example of application is based on data published by [Raywade et al.(2020)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.3551R/abstract). The Authors collected [fast-radio bursts (FRBs)](https://en.wikipedia.org/wiki/Fast_radio_burst) from the event **FRB 121102** in order to identify a possible periodicity in the occurrence these events.


```julia
using CairoMakie
using Downloads: download
using Literate
using ZipStreams
```

## Download the dataset
***
First of all, let's download the table with the events selected in this study. It is available as supplemental material at the site of the journal.

```julia
url = "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mnras/495/4/10.1093_mnras_staa1237/1/staa1237_supplemental_files.zip?Expires=1766572517&Signature=adGzR1A9rTR87jcRkYgNoTttES089LNZfiERXSlFMe7PXiETHLtRVleSsvn8CXCEbRFdcPqTq30FYJkocKHpV11yZ395ny1j4zIvYP5uXRrZ47vqtyNw2lF0uHlQPtGL3PI0yFF0kjLEBNcVpfVs1igA753JRkIocmyrjdodoY03mBOVcvlv44riSBgEjLba5Jv-YgpfCnuq2J2NkDlh~N-GWM8QZu3frI4IhlMbbC69hai8Q9a4hc6OTHh1T7pqYcx-q78JT5tZBRWUUjeUU~cBl44OjZsALEFrYS4wyoQxVqunw9LX8ti-sTNCvn4QCb42cqWtw8KsqiXAikYR9Q__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA"
```

```julia
data = Vector{Float64}()

zipsource(download(url)) do source   # context management of sources with "do" syntax
    for f in source                  # iterate through files in an archive
        println(info(f).name)
        if info(f).name == "All_dets_paper.txt"
            open(info(f).name) do fdt
                l = readlines(fdt)
                for i in l
                    push!(data,parse(Float64,split(i)[2]))
                end
            end
        end
    end
end
```

```julia
println(data)
```



