using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using DelimitedFiles
using HDF5
using Plots
using LaTeXStrings
include("scripts/utils.jl")
include("scripts/write_hdf5.jl")
include("scripts/write_correlatormatrix.jl")
include("scripts/eigenvalues.jl")
include("scripts/massplots.jl")
include("scripts/plotting.jl")
include("scripts/tables.jl")
include("scripts/spectrumplot.jl")
include("scripts/tex_tables.jl")
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=12, tickfontsize=12, labelfontsize=14, titlefontsize=14,  markersize=5)

start_from_logs  = false
write_correlator = false
write_any_hdf5   = false

Nsmear = collect(0:10:80)

logpath        = "/home/fabian/Downloads/DiaLData/measurements"
hdf5path       = "/home/fabian/Downloads/hdf5out_modified"
paramter_path  = "input"
corrfitterpath = "output/fitresults"
tablepath      = "output/tables"
plotpath       = "output/plots"

parameterfile      = joinpath(paramter_path,"parameters_smeared.csv")
parameters_fitting = joinpath(paramter_path,"parameters_corrfitter.csv")
parameters_gevp    = joinpath(paramter_path,"parameters_gevp.csv")
ispath(corrfitterpath) || mkpath(corrfitterpath)

start_from_logs*write_any_hdf5  && main_write_hdf5_logs(Nsmear,logpath,hdf5path,parameterfile)
write_correlator*write_any_hdf5 && main_write_correlator_matrices(Nsmear,hdf5path)
write_any_hdf5 && write_eigenvalues(parameters_gevp,hdf5path)
#all_effective_mass_plots(hdf5path,parameters_gevp)

function run_corrfitter(parameters_fitting,hdf5path;resample=false)
    resample = resample ? "True" : "False"
    args = `$(abspath(parameters_fitting)) $(abspath(hdf5path)) $(abspath(corrfitterpath)) $resample`
    try
        run(`python3 scripts/fitting_eigenvalues.py $args`)
    catch
        run(`python  scripts/fitting_eigenvalues.py $args`)
    end
end

run_corrfitter(parameters_fitting,hdf5path;resample=false)
plot_all_masses_with_fitting(parameters_gevp,parameters_fitting,corrfitterpath,hdf5path,plotpath;only_singlet=true)
write_all_tables(Nsmear,parameters_gevp,parameters_fitting,corrfitterpath,tablepath)
write_tex_tables(tablepath,tablepath)
plot_spectrum(tablepath,plotpath)


