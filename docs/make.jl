using Documenter, GibbsTypePriors

makedocs(
    modules = [GibbsTypePriors],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "konkam",
    sitename = "GibbsTypePriors.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/konkam/GibbsTypePriors.jl.git",
    push_preview = true
)
