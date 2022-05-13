using Documenter, ErrorsInVariablesExtremes

makedocs(sitename="ErrorsInVariablesExtremes.jl",
		 pages = ["index.md", "functions.md"])

deploydocs(repo = "github.com/JuliaExtremes/ErrorsInVariablesExtremes.jl.git")
