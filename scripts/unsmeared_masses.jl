using Pkg; Pkg.activate(".")
using DelimitedFiles
using HiRepParsing
using MixedRepSinglets
using HDF5
using Plots
gr(frame=:box)
plotlyjs(frame=:box)
include("utils.jl")

h5file   = "/home/fabian/Downloads/b55_tests.hdf5"
datapath = "/home/fabian/Dokumente/DataDiaL/"
datapath = "/home/fabian/Documents/Physics/Data/DataDiaL/"

isfile(h5file) && rm(h5file)
ispath("output") || mkpath("output")

path = joinpath(datapath,"measurementsTests")
parameterfile655 = "input/parameters_b6p55.csv"
main_write_hdf5_logs(path,h5file,parameterfile655)

path = joinpath(datapath,"measurements")
parameterfile645 = "input/parameters_b6p45.csv"
main_write_hdf5_logs(path,h5file,parameterfile645)

parameters_fitting = "input/parameters_fitting.csv"
run_corrfitter(parameters_fitting,h5file,"output")

parameters_fitting_fpi = "input/parameters_fitting_fpi.csv"
run_corrfitter(parameters_fitting_fpi,h5file,"output",fpi=true)

#results  = readdlm("output/corrfitter_results.csv",';')
#fitparam = readdlm("input/parameters_fitting.csv",';',skipstart=1)
#_effective_mass_plots(results,fitparam)

results  = readdlm("output/corrfitter_fpi_results.csv",';')
fitparam = readdlm("input/parameters_fitting_fpi.csv",';',skipstart=1)
_effective_mass_plots(results,fitparam)