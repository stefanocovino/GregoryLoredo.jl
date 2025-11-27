```@julia
using CairoMakie
using Downloads: download
using Format
using ZipStreams
```

# A real life example of the Gregory & Loredo (1992) algorithm
***

This example of application of the GL92 algorithm is based on data published by [Raywade et al.(2020)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.3551R/abstract). The Authors collected [fast-radio bursts (FRBs)](https://en.wikipedia.org/wiki/Fast_radio_burst) from the event **FRB 121102** in order to identify a possible periodicity in the occurrence these events.




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

## Visualize the dataset
***

Let's plot an histogram showing the FRB detection epoch distribution in our dataset.

```julia
f = Figure()

ax = Axis(f[1,1],
    xlabel = "Time (MJD)",
    ylabel = "N",
    )

hist!(data)

f
```
![Histogram](histdata.png)


```julia
println("Data consists of ", length(data), " events.")
```

```julia
Data consists of 234 events.
```

We need to define a frequency range for the analysis. We can provide a range defined manually, or rely on an automatic generation baased on data length and range. We can have a linear or logarithmic spacing.


```julia
wr = FrequencyRange(data,logspacing=true)

whi = maximum(wr)
wlo = minimum(wr)

printfmtln("Minimum period: {:.1f} days, maximum period: {:.1f} days, number of periods: {:d}",2π/whi,2π/wlo,length(wr))
```

```julia
Minimum period: 10.7 days, maximum period: 250.3 days, number of periods: 117
```

## Model parameter
***

The
