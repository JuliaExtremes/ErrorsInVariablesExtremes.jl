var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/#Methods","page":"Functions","title":"Methods","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"Modules = [ErrorsInVariablesExtremes]\nPrivate = false\nOrder = [:function]\nPages = [\"src/structures/pseudodata.jl\",\n\t\"src/structures/pseudoensemble.jl\",\n\t\"src/structures/pseudomaximaeva.jl\"\n\t]","category":"page"},{"location":"functions/#Distributions.logpdf-Tuple{Pseudodata, Vector{var\"#s5\"} where var\"#s5\"<:Real}","page":"Functions","title":"Distributions.logpdf","text":"logpdf(pdata::Pseudodata, y::Vector{<:Real})\n\nCompute the log density of each of the potential data y according to the distributions in pdata.\n\n\n\n\n\n","category":"method"},{"location":"functions/#Distributions.pdf-Tuple{Pseudodata, Vector{var\"#s5\"} where var\"#s5\"<:Real}","page":"Functions","title":"Distributions.pdf","text":"pdf(pdata::Pseudodata, y::Vector{<:Real})\n\nCompute the density of each of the potential data y according to the distributions in pdata.\n\n\n\n\n\n","category":"method"},{"location":"functions/#Base.convert-Tuple{Type{DataFrames.DataFrame}, Pseudoensemble}","page":"Functions","title":"Base.convert","text":"convert(::Type{DataFrame}, pensemble::Pseudoensemble)\n\nConvert the Pseudoensemble type to a DataFrame\n\n\n\n\n\n","category":"method"},{"location":"functions/#Distributions.logpdf-Tuple{Pseudoensemble, Vector{var\"#s34\"} where var\"#s34\"<:Real}","page":"Functions","title":"Distributions.logpdf","text":"logpdf(pensemble::Pseudoensemble, y::Vector{<:Real})\n\nCompute the logpdf of each of the potential data y according to the distributions in pensemble.\n\nDetails\n\nIndependance is assumed between the members of the pseudoensemble.value, i.e. the sum of the logpdf of each member is taken.\n\n\n\n\n\n","category":"method"},{"location":"functions/#Distributions.pdf-Tuple{Pseudoensemble, Vector{var\"#s35\"} where var\"#s35\"<:Real}","page":"Functions","title":"Distributions.pdf","text":"pdf(pensemble::Pseudoensemble, y::Vector{<:Real})\n\nCompute the pdf of each of the potential data y according to the distributions in pensemble.\n\nDetails\n\nIndependance is assumed between the members of the pseudoensemble.value, i.e. the product of the pdf of each member is taken.\n\n\n\n\n\n","category":"method"},{"location":"functions/#ErrorsInVariablesExtremes.ensemblemean-Tuple{Pseudoensemble}","page":"Functions","title":"ErrorsInVariablesExtremes.ensemblemean","text":"ensemblemean(pdata::Pseudoensemble)\n\nCompute the ensemble mean for each year.\n\n\n\n\n\n","category":"method"},{"location":"functions/#Extremes.loglike-Tuple{PseudoMaximaEVA, Vector{var\"#s34\"} where var\"#s34\"<:Real, AbstractVector{var\"#s33\"} where var\"#s33\"<:Real}","page":"Functions","title":"Extremes.loglike","text":"function loglike(fm::PseudoMaximaEVA, y::Vector{<:Real}, θ::AbstractVector{<:Real})\n\nCompute the loglikelihood of the ErrorsInVariables extreme value model fm given the maxima y and the GEV parameters θ.\n\n\n\n\n\n","category":"method"},{"location":"functions/#Extremes.loglike-Tuple{PseudoMaximaEVA}","page":"Functions","title":"Extremes.loglike","text":"function loglike(fm::PseudoMaximaEVA)\n\nCompute the loglikelihood of the ErrorsInVariables extreme value model fm for all MCMC iterations.\n\n\n\n\n\n","category":"method"},{"location":"functions/#Mamba.dic-Tuple{PseudoMaximaEVA}","page":"Functions","title":"Mamba.dic","text":"dic(fm::PseudoMaximaEVA)\n\nCompute the Deviance Information Criterion (DIC) described by Gelman et al. (2013) for the PseudoMaximaEVA model fm.\n\nDetails\n\nReference: Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A. & Rubin, D.B. (2013). Bayesian Data Analysis (3rd ed.). Chapman and Hall/CRC. https://doi.org/10.1201/b16018\n\n\n\n\n\n","category":"method"},{"location":"functions/#Types","page":"Functions","title":"Types","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"Modules = [ErrorsInVariablesExtremes]\nPrivate = false\nOrder = [:type]","category":"page"},{"location":"functions/#ErrorsInVariablesExtremes.PseudoMaximaEVA","page":"Functions","title":"ErrorsInVariablesExtremes.PseudoMaximaEVA","text":"PseudoMaximaEVA(pseudodata::Pseudoensemble, model::EVA, sim::Mamba.Chains)\n\nConstruct a PseudoMaximaEVA type.\n\nDetails\n\nEncapsulates the probability distributions for each of the unobserved data for an ensemble of pseudodata, the extreme value model and a sample from the parameters posterior distribution.\n\nTODO: Verify if all member has the same year vector.\n\n\n\n\n\n","category":"type"},{"location":"functions/#ErrorsInVariablesExtremes.Pseudodata","page":"Functions","title":"ErrorsInVariablesExtremes.Pseudodata","text":"Pseudodata(name::String, value::Vector{<:UnivariateDistribution})\n\nConstruct a Pseudodata type.\n\nDetails\n\nEncapsulates the probability distributions for each of the unobserved data. The data are not directly observed but their distributions are known. The distributions can be different for each of the data.\n\n\n\n\n\n","category":"type"},{"location":"functions/#ErrorsInVariablesExtremes.Pseudoensemble","page":"Functions","title":"ErrorsInVariablesExtremes.Pseudoensemble","text":"Pseudoensemble(name::String, value::Vector{Pseudodata})\n\nConstruct a Pseudoensemble type.\n\nDetails\n\nEncapsulates the probability distributions for each of the unobserved data for an ensemble of pseudodata. For the moment, independance of each Pseudodata elements is assumed.\n\nTODO: Verify if all member has the same year vector.\n\n\n\n\n\n","category":"type"},{"location":"tutorial/pseudoensemble/#Pseudoensemble","page":"Pseudoensemble","title":"Pseudoensemble","text":"","category":"section"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"The Pseudoensemble data structure is illustrated using the unobserved riverflows at an ungauged section of the Chaudière River in Québec from 1960 to 2020.","category":"page"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"using ErrorsInVariablesExtremes, Extremes, Dates, DataFrames, Distributions, Gadfly","category":"page"},{"location":"tutorial/pseudoensemble/#Load-the-data","page":"Pseudoensemble","title":"Load the data","text":"","category":"section"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"Loading the pseudo-observations into an object of type Pseudoensemble:","category":"page"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"filename = \"../../../test/data/A2020_Analyse_Historique_QMA_SLSO00003.nc\"\npensemble = ErrorsInVariablesExtremes.load_discharge_distribution(filename)","category":"page"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"The vector of years can be extracted:","category":"page"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"years = pensemble.value[1].year;","category":"page"},{"location":"tutorial/pseudoensemble/#Convert-to-DataFrame","page":"Pseudoensemble","title":"Convert to DataFrame","text":"","category":"section"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"Converting Pseudoensemble to DataFrame:","category":"page"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"df = convert(DataFrame, pensemble)\nfirst(df,5)","category":"page"},{"location":"tutorial/pseudoensemble/#Compute-the-mean-and-the-1-α-confidence-interval-for-each-of-the-configurations","page":"Pseudoensemble","title":"Compute the mean and the 1-α confidence interval for each of the configurations","text":"","category":"section"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"α = 0.5\n\ndf[:, :y] = mean.(df.Distribution)\ndf[:, :ymin] = quantile.(df.Distribution, α/2)\ndf[:, :ymax] = quantile.(df.Distribution, 1-α/2)\n\nfirst(df, 5)","category":"page"},{"location":"tutorial/pseudoensemble/#Display-the-discharge-of-all-configurations","page":"Pseudoensemble","title":"Display the discharge of all configurations","text":"","category":"section"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"set_default_plot_size(12cm, 8cm)\nplot(df, x=:Year, y=:y, color=:Configuration, Geom.line,\n    ymin=:ymin, ymax=:ymax, Geom.ribbon,\n    Guide.ylabel(\"Discharge [m³/s]\"))","category":"page"},{"location":"tutorial/pseudoensemble/#Display-the-discharge-of-one-specific-configuration","page":"Pseudoensemble","title":"Display the discharge of one specific configuration","text":"","category":"section"},{"location":"tutorial/pseudoensemble/","page":"Pseudoensemble","title":"Pseudoensemble","text":"config = \"MG24HQ\"\n\ndf2 = filter(row -> row.Configuration ==config, df)\n\nset_default_plot_size(12cm, 8cm)\nplot(df2, x=:Year, y=:y, Geom.line,\n    ymin=:ymin, ymax=:ymax, Geom.ribbon,\n    Guide.ylabel(\"Discharge [m³/s]\"))","category":"page"},{"location":"tutorial/stationary_model/#Stationary","page":"Stationary extremes model","title":"Stationary","text":"","category":"section"},{"location":"tutorial/stationary_model/","page":"Stationary extremes model","title":"Stationary extremes model","text":"using ErrorsInVariablesExtremes, Extremes, Dates, DataFrames, Distributions, Gadfly, Cairo, Fontconfig","category":"page"},{"location":"tutorial/stationary_model/#Load-the-data","page":"Stationary extremes model","title":"Load the data","text":"","category":"section"},{"location":"tutorial/stationary_model/","page":"Stationary extremes model","title":"Stationary extremes model","text":"Loading the pseudo-observations into an object of type Pseudoensemble:","category":"page"},{"location":"tutorial/stationary_model/","page":"Stationary extremes model","title":"Stationary extremes model","text":"filename = \"../../../test/data/A2020_Analyse_Historique_QMA_SLSO00003.nc\"\npensemble = ErrorsInVariablesExtremes.load_discharge_distribution(filename)","category":"page"},{"location":"tutorial/stationary_model/","page":"Stationary extremes model","title":"Stationary extremes model","text":"The vector of years can be extracted:","category":"page"},{"location":"tutorial/stationary_model/","page":"Stationary extremes model","title":"Stationary extremes model","text":"years = pensemble.value[1].year;","category":"page"},{"location":"tutorial/stationary_model/#MCMC","page":"Stationary extremes model","title":"MCMC","text":"","category":"section"},{"location":"tutorial/stationary_model/","page":"Stationary extremes model","title":"Stationary extremes model","text":"res = gevfitbayes(pensemble)","category":"page"},{"location":"tutorial/stationary_model/#Display-of-annual-maximum-estimates","page":"Stationary extremes model","title":"Display of annual maximum estimates","text":"","category":"section"},{"location":"tutorial/stationary_model/","page":"Stationary extremes model","title":"Stationary extremes model","text":"\n# Compute the estimates\nŷ = vec(mean( res.maxima[:, :, 1].value, dims=1))\n\n# Compute the 95% credible interval\nymin = [quantile(vec(res.maxima[:, j, 1].value), .025) for j in 1:length(ŷ)]\nymax= [quantile(vec(res.maxima[:, j, 1].value), .975) for j in 1:length(ŷ)]\n\ndf = DataFrame(Year = years, Discharge=ŷ, ymin = ymin, ymax=ymax)\n\n# plot the ponctual estimates and the intervals\nset_default_plot_size(12cm, 8cm)\nplot(df, x=:Year, y=:Discharge, Geom.line, Geom.point,\n    ymin=ymin, ymax=ymax, Geom.ribbon\n)","category":"page"},{"location":"tutorial/stationary_model/#Fit-of-the-GEV-law-adjusted-to-the-estimates-of-the-annual-maxima","page":"Stationary extremes model","title":"Fit of the GEV law adjusted to the estimates of the annual maxima","text":"","category":"section"},{"location":"tutorial/stationary_model/","page":"Stationary extremes model","title":"Stationary extremes model","text":"set_default_plot_size(21cm ,16cm)\ndiagnosticplots(res.parameters)","category":"page"},{"location":"tutorial/#Tutorial","page":"Getting started","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Getting started","title":"Getting started","text":"This tutorial shows the functionalities of ErrorsInVariablesExtremes.jl. They are illustrated by producing an Extreme Value Analysis of unobserved riverflows at ungauged section of the Chaudière River.","category":"page"},{"location":"tutorial/","page":"Getting started","title":"Getting started","text":"Before executing this tutorial, make sure to have installed the following packages:","category":"page"},{"location":"tutorial/","page":"Getting started","title":"Getting started","text":"ErrorsInVariablesExtremes.jl (of course)\nExtremes.jl (for more Extreme Value Analysis tools)\nDataFrames.jl (for using the DataFrame type)\nDistributions.jl (for using probability distribution objects)\nGadfly.jl (for plotting)\nCSV.jl (to load CSV files)","category":"page"},{"location":"tutorial/","page":"Getting started","title":"Getting started","text":"and import them using the following command:","category":"page"},{"location":"tutorial/","page":"Getting started","title":"Getting started","text":"using ErrorsInVariablesExtremes, Extremes, Dates, DataFrames, Distributions, Gadfly, CSV","category":"page"},{"location":"#Extreme-value-analysis-for-error-in-variables-models.","page":"Extreme value analysis for error in variables models.","title":"Extreme value analysis for error in variables models.","text":"","category":"section"},{"location":"tutorial/non_stationary_model/#Stationary","page":"Non-stationary extremes model","title":"Stationary","text":"","category":"section"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"using ErrorsInVariablesExtremes, Extremes, Dates, DataFrames, Distributions, Gadfly, Cairo, Fontconfig, CSV","category":"page"},{"location":"tutorial/non_stationary_model/#Load-the-data","page":"Non-stationary extremes model","title":"Load the data","text":"","category":"section"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"Loading the pseudo-observations into an object of type Pseudoensemble:","category":"page"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"filename = \"../../../test/data/A2020_Analyse_Historique_QMA_SLSO00003.nc\"\npensemble = ErrorsInVariablesExtremes.load_discharge_distribution(filename)","category":"page"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"The vector of years can be extracted:","category":"page"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"years = pensemble.value[1].year;","category":"page"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"Loading the covariable:","category":"page"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"co2data = CSV.read(\"../../../test/data/RCPdata.csv\", DataFrame);\nfilter!(row->row.Year in pensemble.value[1].year, co2data)\nlocationcov= Extremes.buildVariables(co2data, [:RCP85])","category":"page"},{"location":"tutorial/non_stationary_model/#MCMC","page":"Non-stationary extremes model","title":"MCMC","text":"","category":"section"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"res = gevfitbayes(pensemble, locationcov = locationcov)","category":"page"},{"location":"tutorial/non_stationary_model/#Display-of-annual-maximum-estimates","page":"Non-stationary extremes model","title":"Display of annual maximum estimates","text":"","category":"section"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"# Compute the estimates\nŷ = vec(mean( res.maxima[:, :, 1].value, dims=1))\n\n# Compute the 95% credible interval\nymin = [quantile(vec(res.maxima[:, j, 1].value), .025) for j in 1:length(ŷ)]\nymax= [quantile(vec(res.maxima[:, j, 1].value), .975) for j in 1:length(ŷ)]\n\ndf = DataFrame(Year = years, Discharge=ŷ, ymin = ymin, ymax=ymax)\n\n# plot the ponctual estimates and the intervals\nset_default_plot_size(12cm, 8cm)\nplot(df, x=:Year, y=:Discharge, Geom.line, Geom.point,\n    ymin=ymin, ymax=ymax, Geom.ribbon\n)","category":"page"},{"location":"tutorial/non_stationary_model/#Fit-of-the-GEV-law-adjusted-to-the-estimates-of-the-annual-maxima","page":"Non-stationary extremes model","title":"Fit of the GEV law adjusted to the estimates of the annual maxima","text":"","category":"section"},{"location":"tutorial/non_stationary_model/","page":"Non-stationary extremes model","title":"Non-stationary extremes model","text":"set_default_plot_size(21cm ,16cm)\ndiagnosticplots(res.parameters)","category":"page"}]
}
