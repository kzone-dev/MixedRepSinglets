using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
using LaTeXStrings
using DelimitedFiles
include("utils.jl")
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=11, tickfontsize=10, labelfontsize=14, markersize=5)

h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M1234_with_resamples_more_bins_v2.hdf5"
parameters = readdlm("input/parameters_gevp.csv",';';skipstart=1)
nrows = first(size(parameters))

for row in 1:nrows

    row > 4 && continue

    ensemble, channel, t0, binsize, deriv, ops = parameters[row,:]
    
    nops = parse.(Int,split(replace(ops,r"[()]"=>""),','))

    β   = h5read(h5eigenvals,joinpath(ensemble,channel,"beta"))
    T,L = h5read(h5eigenvals,joinpath(ensemble,channel,"lattice"))[1:2]
    mf  = h5read(h5eigenvals,joinpath(ensemble,channel,"quarkmasses_fundamental"))[1]
    mas = h5read(h5eigenvals,joinpath(ensemble,channel,"quarkmasses_antisymmetric"))[1]
    meff     = h5read(h5eigenvals,joinpath(ensemble,channel,"meff"))
    eigvals  = h5read(h5eigenvals,joinpath(ensemble,channel,"eigvals"))
    Δmeff    = h5read(h5eigenvals,joinpath(ensemble,channel,"Delta_meff"))
    Δeigvals = h5read(h5eigenvals,joinpath(ensemble,channel,"Delta_eigvals"))

    plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=2)
    plot!(plt1, ylims=(0.3,1.2))
    display(plt1)
end

