using Documenter, ErrorsInVariablesExtremes

makedocs(sitename="ErrorsInVariablesExtremes.jl",
		 pages = ["index.md",
				  "Tutorial" =>["Getting started" => "tutorial/index.md",
								"Pseudoensemble" => "tutorial/pseudoensemble.md",
								"Stationary extremes model" => "tutorial/stationary_model.md",
								"Non-stationary extremes model" => "tutorial/non_stationary_model.md"],
				  "functions.md"])

deploydocs(repo = "github.com/JuliaExtremes/ErrorsInVariablesExtremes.jl.git",
		   branch = "gh-pages",
    	   devbranch = "dev",
    	   devurl = "dev")
