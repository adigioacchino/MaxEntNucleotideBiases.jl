using MaxEntNucleotideBiases
using Documenter

DocMeta.setdocmeta!(MaxEntNucleotideBiases, :DocTestSetup, :(using MaxEntNucleotideBiases); recursive=true)

makedocs(;
    modules=[MaxEntNucleotideBiases],
    authors="Andrea Di Gioacchino",
    repo="https://github.com/adigioacchino/MaxEntNucleotideBiases.jl/blob/{commit}{path}#{line}",
    sitename="MaxEntNucleotideBiases.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://adigioacchino.github.io/MaxEntNucleotideBiases.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/adigioacchino/MaxEntNucleotideBiases.jl",
    devbranch="main",
)
