
using CSV, DataFrames, Dates, NetCDF
using Distributions, Extremes, LinearAlgebra, Mamba, Random, Statistics
using Gadfly
using ProgressMeter

using Test

using ErrorsInVariablesExtremes


filename = "/Users/jalbert/Dropbox (MAGI)/Files/Research/2021/TransfertInfoCrue/data/A2020_Analyse_Historique_QMA_SLSO00003.nc"
pensemble = load_discharge_distribution(filename);

res = gevfitbayes(pensemble)