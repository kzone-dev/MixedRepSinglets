using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using LatticeUtils
using DelimitedFiles
using HDF5
using Plots
using LaTeXStrings
using Statistics
using LinearAlgebra
include("scripts/utils.jl")
include("scripts/write_hdf5.jl")
include("scripts/plotting.jl")
pgfplotsx(legend=:topright, frame=:box, legendfontsize=14, tickfontsize=14, labelfontsize=14, titlefontsize=16,  markersize=5)

parameterfile      = joinpath(paramter_path,"parameters_smeared.csv")

ispath(hdf5file_path) || mkpath(hdf5file_path) 

NsmearFUN = collect(0:50:200)
NsmearAS  = collect(0:60:180)

channels=["g5","id"]
write_correlator   && main_write_correlator_matrices(NsmearFUN,NsmearAS,hdf5file_path, ensemble)