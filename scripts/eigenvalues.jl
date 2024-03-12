using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HDF5
using DelimitedFiles
include("utils.jl")

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M1234.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M1234_with_resamples_more_bins_v2.hdf5"

parameters = readdlm("input/parameters_gevp.csv",';';skipstart=1)

for row in eachrow(parameters)

    ensemble, channel, t0, binsize, deriv, ops = row
    nops = parse.(Int,split(replace(ops,r"[()]"=>""),','))
    
    matrixname ="correlation_matrix_$channel"
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,matrixname))
    correlation_matrix = correlation_matrix[nops,nops,:,:]
    
    write_eigenvalues_and_effective_masses(correlation_matrix,h5eigenvals,h5corrs,ensemble,channel;t0,binsize,deriv,resamples=true)
end




