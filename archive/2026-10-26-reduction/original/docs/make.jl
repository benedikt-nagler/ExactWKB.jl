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
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "The cubic, end to end" => "tutorial_cubic.md",
            "The double well" => "tutorial_double_well.md",
            "Seiberg-Witten SU(2)" => "tutorial_sw.md",
        ],
        "Manual" => [
            "Problems and turning points" => "problems.md",
            "The WKB expansion" => "wkb.md",
            "Stokes graphs" => "stokes.md",
            "Spectra" => "spectra.md",
            "Wall-crossing" => "ddp.md",
            "The cluster bridge" => "bridge.md",
            "Thermodynamic Bethe ansatz" => "tba.md",
            "Seiberg-Witten geometry" => "seiberg_witten.md",
            "Plotting" => "plotting.md",
        ],
        "Errors" => "errors.md",
        "API index" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/benedikt-nagler/ExactWKB.jl", devbranch = "main")
