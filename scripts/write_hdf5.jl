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

path = "/home/fabian/Dokumente/DataDiaL/measurementsTests"
h5file = "/home/fabian/Dokumente/Physics/Analysis/HiRepHadrons/b55_tests.hdf5"
parameterfile = "input/parameters.csv"
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
