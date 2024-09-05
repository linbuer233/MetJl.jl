using MetJl
using Documenter

DocMeta.setdocmeta!(MetJl, :DocTestSetup, :(using MetJl); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

makedocs(;
    modules = [MetJl],
    authors = "linziyang <linziyang233@163.com> and contributors",
    repo = "https://github.com/linbuer233/MetJl.jl/blob/{commit}{path}#{line}",
    sitename = "MetJl.jl",
    format = Documenter.HTML(; canonical = "https://linbuer233.github.io/MetJl.jl"),
    pages = [
        "index.md"
        [
            file for file in readdir(joinpath(@__DIR__, "src")) if
            file != "index.md" && splitext(file)[2] == ".md"
        ]
    ],
)

deploydocs(; repo = "github.com/linbuer233/MetJl.jl")
