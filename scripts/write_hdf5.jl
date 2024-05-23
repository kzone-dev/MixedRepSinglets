using DelimitedFiles
using HiRepParsing

function main_write_hdf5_logs(path,h5file,parameterfile)
    input = readdlm(parameterfile,';',skipstart=1)
    for prm in eachrow(input)
        
        dir, file, typeCONN, rep, name = prm
        fileCONN = joinpath(path,dir,"out",file)

        @show fileCONN   
        writehdf5_spectrum(fileCONN,h5file,typeCONN,h5group="$name/$rep/CONN")
    end    
end

h5file   = "/home/fabian/Downloads/b55_tests.hdf5"
datapath = "/home/fabian/Dokumente/DataDiaL/"
datapath = "/home/fabian/Documents/DataDiaL/"

isfile(h5file) && rm(h5file)
ispath("output") || mkpath("output")

path = joinpath(datapath,"measurementsTests")
parameterfile = "input/parameters_b6p55.csv"
main_write_hdf5_logs(path,h5file,parameterfile)

path = joinpath(datapath,"measurements")
parameterfile = "input/parameters_b6p45.csv"
main_write_hdf5_logs(path,h5file,parameterfile)

function run_corrfitter(parameters,hdf5file,outdir)
    args = `$(abspath(parameters)) $(abspath(hdf5file)) $(abspath(outdir))`
    try
        run(`python3 $(abspath("scripts/fitting_eigenvalues.py")) $args`)
    catch
        run(`python  $(abspath("scripts/fitting_eigenvalues.py")) $args`)
    end
end

parameters_fitting = "input/parameters_fitting.csv"
run_corrfitter(parameters_fitting,h5file,"output")
