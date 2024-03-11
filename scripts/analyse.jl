using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
using LaTeXStrings
using DelimitedFiles
include("utils.jl")
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=11, tickfontsize=10, labelfontsize=14, markersize=5)

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M1234.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M1234.hdf5"

parameters = readdlm("input/parameters_gevp.csv",';';skipstart=1)

for row in eachrow(parameters)

    ensemble, channel, t0, binsize, deriv, ops = row
    nops = parse.(Int,split(replace(ops,r"[()]"=>""),','))
    
    matrixname ="correlation_matrix_$channel"
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,matrixname))
    correlation_matrix = correlation_matrix[nops,nops,:,:]
    
    write_eigenvalues_and_effective_masses(correlation_matrix,h5eigenvals,h5corrs,ensemble,channel;t0,binsize,deriv)
end




