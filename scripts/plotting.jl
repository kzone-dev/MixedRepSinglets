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
parameters_fitting = readdlm("input/parameters_corrfitter.csv",';';skipstart=1)
corrfitter_results = readdlm("output/corrfitter_results.csv",';';skipstart=0)

#check that the number of datasets match
@assert first(size(parameters)) == first(size(parameters_fitting)) == first(size(corrfitter_results)) 
nrows = first(size(parameters))

function channel_tags(channel)
    isequal(channel,"g5_singlet")        && return L"$J^P = 0^-$(singlet)"
    isequal(channel,"g5_nonsinglet_FUN") && return L"$J^P = 0^-$(f)"
    isequal(channel,"g5_nonsinglet_AS")  && return L"$J^P = 0^-$(as)"
    isequal(channel,"g1_nonsinglet_FUN") && return L"$J^P = 1^+$(f)"
    isequal(channel,"g1_nonsinglet_AS")  && return L"$J^P = 1^+$(as)"
    return ""
end

for row in 1:nrows

    row > 4 && continue

    ensemble, channel, t0, binsize, deriv, ops = parameters[row,:]
    ensemble, channel, tmin1, tmin2, tmax1, tmax2, tp, Nmax  = parameters_fitting[row,:]
    ensemble, channel, T, L, mf, mas, beta, E0, ΔE0, E1, ΔE1, χ2dof0, χ2dof0  = corrfitter_results[row,:]

    nops = parse.(Int,split(replace(ops,r"[()]"=>""),','))

    β   = h5read(h5eigenvals,joinpath(ensemble,channel,"beta"))
    T,L = h5read(h5eigenvals,joinpath(ensemble,channel,"lattice"))[1:2]
    mf  = h5read(h5eigenvals,joinpath(ensemble,channel,"quarkmasses_fundamental"))[1]
    mas = h5read(h5eigenvals,joinpath(ensemble,channel,"quarkmasses_antisymmetric"))[1]
    meff     = h5read(h5eigenvals,joinpath(ensemble,channel,"meff"))
    eigvals  = h5read(h5eigenvals,joinpath(ensemble,channel,"eigvals"))
    Δmeff    = h5read(h5eigenvals,joinpath(ensemble,channel,"Delta_meff"))
    Δeigvals = h5read(h5eigenvals,joinpath(ensemble,channel,"Delta_eigvals"))

    plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=2,tag=channel_tags(channel))
    add_fit_range!(plt1,tmin1,tmax1,E0,ΔE0;label="")
    add_fit_range!(plt1,tmin2,tmax2,E1,ΔE1;label="")
    
    @show E1,ΔE1

    plot!(plt1, ylims=(0.8*E0,1.2*E0))
    plot!(plt1, ylims=(0.3,1.2))
    display(plt1)
end

