using Documenter, ErrorsInVariablesExtremes, Extremes, DataFrames

makedocs(sitename="ErrorsInVariablesExtremes.jl",
		 pages = ["index.md",
				  "Tutorial" =>["Getting started" => "tutorial/index.md",
								"Unobserved riverflows analysis" => "tutorial/riverflows.md"],
				  "functions.md"])

deploydocs(repo = "github.com/JuliaExtremes/ErrorsInVariablesExtremes.jl.git",
		   branch = "gh-pages",
    	   devbranch = "dev",
    	   devurl = "dev")
