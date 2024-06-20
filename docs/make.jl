using Documenter, PyramidScheme

makedocs(
    modules = [PyramidScheme],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Felix Cremer, Fabian Gans",
    sitename = "PyramidScheme.jl",
    pages = Any[
            "index.md",
            "interface.md",
    ]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/JuliaDataCubes/PyramidScheme.jl.git",
    push_preview = true
)
