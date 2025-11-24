# From [Rajwade et al.(2020)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.3551R/abstract) we downloaded a file with the arrival time, expressed in MJD, of FRBs detected from **FRB 121102**. For any detail about data collection and processing the quoted paper is the best reference.

using CairoMakie
using Downloads: download
using Literate
using ZipStreams

# ## Download the dataset
# ***
# Firs of all, let's download the table with the events selected in this study. It is available as supplemental material at the site of the journal.

url = "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mnras/495/4/10.1093_mnras_staa1237/1/staa1237_supplemental_files.zip?Expires=1766572517&Signature=adGzR1A9rTR87jcRkYgNoTttES089LNZfiERXSlFMe7PXiETHLtRVleSsvn8CXCEbRFdcPqTq30FYJkocKHpV11yZ395ny1j4zIvYP5uXRrZ47vqtyNw2lF0uHlQPtGL3PI0yFF0kjLEBNcVpfVs1igA753JRkIocmyrjdodoY03mBOVcvlv44riSBgEjLba5Jv-YgpfCnuq2J2NkDlh~N-GWM8QZu3frI4IhlMbbC69hai8Q9a4hc6OTHh1T7pqYcx-q78JT5tZBRWUUjeUU~cBl44OjZsALEFrYS4wyoQxVqunw9LX8ti-sTNCvn4QCb42cqWtw8KsqiXAikYR9Q__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA"


println(download(url))

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









#Literate.markdown("examples/README.jl", "."; flavor = Literate.CommonMarkFlavor())
#Literate.notebook(inputfile, outputdir=pwd(); config::AbstractDict=Dict(), kwargs...)
