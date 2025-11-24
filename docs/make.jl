using GregoryLoredo
using Documenter
#using Literate



# generate examples
#EXAMPLE = "src/Example.jl"
#OUTPUTJ = "build/Example.ipynb"
#OUTPUTMD = "build/Example.md"

#Literate.markdown(EXAMPLE, OUTPUTMD)
#Literate.notebook(EXAMPLE, OUTPUTJ)



DocMeta.setdocmeta!(GregoryLoredo, :DocTestSetup, :(using GregoryLoredo); recursive=true)

makedocs(;
    modules=[GregoryLoredo],
    authors="Stefano Covino <stefano.covino@inaf.it> and contributors",
    sitename="GregoryLoredo.jl",
    format=Documenter.HTML(;
        canonical="https://stefanocovino.github.io/GregoryLoredo.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        #"Example" => "Example.jl"
    ],
)

deploydocs(;
    repo="github.com/stefanocovino/GregoryLoredo.jl",
    devbranch="main",
)
