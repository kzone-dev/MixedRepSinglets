using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using DelimitedFiles
using HDF5
using Plots
using LaTeXStrings
using Statistics
using LinearAlgebra
using LsqFit
include("scripts/utils.jl")
include("scripts/write_hdf5.jl")
include("scripts/write_correlatormatrix.jl")
include("scripts/eigenvalues.jl")
include("scripts/massplots.jl")
include("scripts/plotting.jl")
include("scripts/tables.jl")
include("scripts/spectrumplot.jl")
include("scripts/tex_tables.jl")
include("scripts/mixing_angle.jl")
pgfplotsx(legend=:topright, frame=:box, legendfontsize=14, tickfontsize=14, labelfontsize=14, titlefontsize=16,  markersize=5)
Nsmear = collect(0:10:80)

corrfitterpath = joinpath(output_path,"fitresults")
tablepath      = joinpath(output_path,"tables")
plotpath       = joinpath(output_path,"plots")

parameterfile      = joinpath(paramter_path,"parameters_smeared.csv")
parameters_fitting = joinpath(paramter_path,"parameters_corrfitter.csv")
parameters_gevp    = joinpath(paramter_path,"parameters_gevp.csv")
gradient_flow_results = joinpath(paramter_path,"gradient_flow_results.csv")

ispath(corrfitterpath) || mkpath(corrfitterpath)
ispath(hdf5file_path) || mkpath(hdf5file_path) 

start_from_logs    && main_write_hdf5_logs(Nsmear,logfiles_path,hdf5file_path,parameterfile;filter_channels=!write_all_channes_to_hdf5)
write_correlator   && main_write_correlator_matrices(Nsmear,hdf5file_path)
write_gevp_results && write_eigenvalues(parameters_gevp,hdf5file_path)

function run_corrfitter(parameters_fitting,hdf5file_path;resample)
    resample = resample ? "True" : "False"
    args = `$(abspath(parameters_fitting)) $(abspath(hdf5file_path)) $(abspath(corrfitterpath)) $resample`
    try
        run(`python3 scripts/fitting_eigenvalues.py $args`)
    catch
        run(`python  scripts/fitting_eigenvalues.py $args`)
    end
end

run_corrfitter(parameters_fitting,hdf5file_path;resample=true)
plot_all_masses_with_fitting(parameters_gevp,parameters_fitting,corrfitterpath,hdf5file_path,plotpath;only_singlet=false)
write_all_tables(Nsmear,parameters_gevp,parameters_fitting,corrfitterpath,tablepath)
write_tex_tables(tablepath,tablepath)
plot_spectrum(tablepath,plotpath,gradient_flow_results)
plot_and_write_mixing_angles(parameters_gevp,hdf5file_path,tablepath,tablepath,plotpath)