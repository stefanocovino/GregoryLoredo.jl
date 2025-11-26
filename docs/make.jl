using GregoryLoredo
using Documenter



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
        "Example" => "Example.md"
    ],
)

deploydocs(;
    repo="github.com/stefanocovino/GregoryLoredo.jl",
    devbranch="main",
)
