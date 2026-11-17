using Documenter
using ExactWKB

DocMeta.setdocmeta!(ExactWKB, :DocTestSetup, :(using ExactWKB); recursive = true)

makedocs(;
    modules = [ExactWKB],
    authors = "Benedikt Nagler <benedikt.nagler@protonmail.com>",
    sitename = "ExactWKB.jl",
    format = Documenter.HTML(;
        canonical = "https://benedikt-nagler.github.io/ExactWKB.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
    ),
    pages = ["Home" => "index.md", "Expansion and period API" => "api.md", "Graphs and numerics" => "numerics.md"],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/benedikt-nagler/ExactWKB.jl", devbranch = "main")
