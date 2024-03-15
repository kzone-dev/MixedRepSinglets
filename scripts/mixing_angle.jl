using MixedRepSinglets
using LaTeXStrings
using DelimitedFiles
using HDF5
using Statistics
using MixedRepSinglets
using Plots
include("utils.jl")
include("eigenvalues.jl")

asin_deriv(x) = +1/sqrt(1-x^2)
acos_deriv(x) = -1/sqrt(1-x^2)
asin_error(x,Δx) = abs(asin_deriv(x))*Δx
acos_error(x,Δx) = abs(acos_deriv(x))*Δx
function plot_effective_mixing_angle(evecs, Δevecs)
    plt = plot(legend=:outerright)
    scatter!(plt, evecs[1,1,:],  ms = 8, yerr = Δevecs[1,1,:],label="(1,1)")
    scatter!(plt, evecs[1,2,:],  ms = 8, yerr = Δevecs[1,2,:],label="(1,2)")
    scatter!(plt, evecs[2,1,:],  ms = 5, yerr = Δevecs[2,1,:],label="(2,1)")
    scatter!(plt, evecs[2,2,:],  ms = 5, yerr = Δevecs[2,2,:],label="(2,2)")
    return plt
end

parameters_gevp = joinpath("input","parameters_gevp.csv")
hdf5path    = "/home/fabian/Downloads/hdf5out_modified"
h5corrs     = joinpath(hdf5path,"singlets_smeared_correlators.hdf5")
h5eigenvals = joinpath(hdf5path,"singlets_smeared_eigenvalues.hdf5")

parameters = readdlm(parameters_gevp,';';skipstart=1)
for row in eachrow(parameters)

    ensemble, channel, t0, binsize, = row[1:4]
    nops, deriv  = [1,10], false
    t0  = 5

    channel == "g5_singlet" || continue    
    matrixname ="correlation_matrix_g5_singlet"
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,matrixname))
    correlation_matrix = correlation_matrix[nops,nops,:,:]
    
    evals, Δevals, meff, Δmeff, evals_jk, evecs, Δevecs, evecs_jk = eigenvalues_eigenvectors_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    plt = plot_effective_mixing_angle(evecs, Δevecs)
    plot!(ylims=(-1.02,1.02))
    display(plt)
end
