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
        canonical="https://eXplorationLogisticsControl.github.io/Sequoia.jl/dev/",
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

# Deployment is done by .github/workflows/Documentation.yml, which publishes docs/build
# to GitHub Pages under dev/.
