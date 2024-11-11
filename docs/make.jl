using Documenter, ErrorsInVariablesExtremes

makedocs(sitename="ErrorsInVariablesExtremes.jl",
		 pages = ["index.md",
				  "Tutorial" =>["Getting started" => "tutorial/index.md",
								"Pseudoensemble" => "tutorial/pseudoensemble.md",
								"Analysis of estimated discharges" => "tutorial/analysis.md"],
				  "functions.md"])

deploydocs(repo = "github.com/JuliaExtremes/ErrorsInVariablesExtremes.jl.git",
		   branch = "gh-pages",
    	   devbranch = "dev",
    	   devurl = "dev")
