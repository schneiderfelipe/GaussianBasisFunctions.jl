using GaussianBasisFunctions
using Documenter

makedocs(;
    modules=[GaussianBasisFunctions],
    authors="Felipe S. S. Schneider <schneider.felipe@posgrad.ufsc.br> and contributors",
    repo="https://github.com/schneiderfelipe/GaussianBasisFunctions.jl/blob/{commit}{path}#L{line}",
    sitename="GaussianBasisFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://schneiderfelipe.github.io/GaussianBasisFunctions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/schneiderfelipe/GaussianBasisFunctions.jl",
)
