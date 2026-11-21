using Documenter
using Sequoia
using JuMP

DocMeta.setdocmeta!(Sequoia, :DocTestSetup, :(using Sequoia); recursive=true)

makedocs(;
    modules=[Sequoia],
    sitename="Sequoia.jl",
    authors="Yuri Shimane",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://eXplorationLogisticsControl.github.io/Sequoia.jl",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "Guide" => [
            "Solving JuMP models" => "guide/jump.md",
            "Custom node oracles" => "guide/oracle.md",
            "Algorithm" => "guide/algorithm.md",
        ],
        "API reference" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/eXplorationLogisticsControl/Sequoia.jl",
    devbranch="main",
    push_preview=true,
)
