using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using DelimitedFiles
using HDF5
using Plots
using LaTeXStrings
include("scripts/utils.jl")
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=10, tickfontsize=10, labelfontsize=12, markersize=5)

start_from_logs  = false
write_correlator = false
write_any_hdf5   = false

Nsmear = collect(0:10:80)

# TODO: Don't expose files, expose paths
logpath        = "/home/fabian/Downloads/DiaLData/measurements"
paramter_path  = "input"
corrfitterpath = "output"
hdf5path       = "output/hdf5files"
tablepath      = "output/tables"
textablepath   = "output/tables"
plotpath       = "plots"

parameterfile      = "input/parameters_smeared.csv"
parameters_fitting = "input/parameters_corrfitter.csv"
parameters_gevp    = "input/parameters_gevp.csv"

h5logfiles    = "output/hdf5files/singlets_smeared.hdf5"
h5corrs       = "output/hdf5files/singlets_smeared_correlators.hdf5"
h5eigenvals   = "output/hdf5files/singlets_smeared_eigenvalues.hdf5"

results_corrfitter    = "output/corrfitter_results.csv"
results_corrfitter_HR = "output/corrfitter_results_HR.csv"

if start_from_logs*write_any_hdf5
    include("scripts/write_hdf5.jl")
    main_write_hdf5_logs(Nsmear,logpath,h5logfiles,parameterfile)
end

if write_correlator*write_any_hdf5
    include("scripts/write_correlatormatrix.jl")
    main_write_correlator_matrices(Nsmear,h5logfiles,h5corrs)
end

if write_any_hdf5
    include("scripts/eigenvalues.jl")
    write_eigenvalues(parameters_gevp,h5corrs,h5eigenvals)

    include("scripts/massplots.jl")
    all_effective_mass_plots(h5eigenvals,parameters_gevp)
end

function run_corrfitter(parameters_fitting,h5eigenvals;resample=false)
    resample = resample ? "True" : "False"
    try
        run(`python3 scripts/fitting_eigenvalues.py $(abspath(parameters_fitting)) $h5eigenvals $resample`)
    catch
        run(`python  scripts/fitting_eigenvalues.py $(abspath(parameters_fitting)) $h5eigenvals $resample`)
    end
end
run_corrfitter(parameters_fitting,h5eigenvals;resample=false)

include("scripts/plotting.jl")
plot_all_masses_with_fitting(parameters_gevp,parameters_fitting,results_corrfitter,h5eigenvals)

include("scripts/tables.jl")
write_all_tables(Nsmear,parameters_gevp,parameters_fitting,results_corrfitter,results_corrfitter_HR)

include("scripts/tex_tables.jl")
write_tex_tables()

include("scripts/spectrumplot.jl")
plot_spectrum("output/tables/table_results_MR.csv")


