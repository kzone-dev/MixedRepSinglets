using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using DelimitedFiles
using HDF5
using Plots
using LaTeXStrings
include("utils.jl")
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=10, tickfontsize=10, labelfontsize=12, markersize=5)

start_from_logs  = false
write_correlator = false

Nsmear        = collect(0:10:80)
logpath       = "/home/fabian/Downloads/DiaLData/measurements"
parameterfile = "input/parameters_smeared.csv"
h5logfiles    = "/home/fabian/Downloads/smeared_singlets_M1234.hdf5"
h5corrs       = "/home/fabian/Downloads/smeared_singlet_correlators_M1234.hdf5"
h5eigenvals   = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M1234_with_resamples_more_bins_v2.hdf5"
gevp_parameterfile = "input/parameters_gevp.csv"

if start_from_logs
    include("write_hdf5.jl")
    main_write_hdf5_logs(Nsmear,logpath,h5logfiles,parameterfile)
end
if write_correlator
    include("write_correlatormatrix.jl")
    main_write_correlator_matrices(Nsmear,h5logfiles,h5corrs)
end
#include("eigenvalues.jl")
#write_eigenvalues(gevp_parameterfile,h5corrs,h5eigenvals)
include("massplots.jl")
all_effective_mass_plots(h5eigenvals,gevp_parameterfile)