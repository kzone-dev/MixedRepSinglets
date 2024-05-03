using DelimitedFiles
using HiRepParsing

function main_write_hdf5_logs(path,hdf5path,parameterfile)
    h5file = joinpath(hdf5path,"b55_tests.hdf5")
    input = readdlm(parameterfile,';',skipstart=1)
    for prm in eachrow(input)

        dir, file, typeCONN, rep, name = prm
        fileCONN = joinpath(path,dir,"out",file)
        
        @show fileCONN   
        writehdf5_spectrum(fileCONN,h5file,typeCONN,h5group="$name/$rep/CONN")
    end    
end

path = "/home/fabian/Dokumente/DataDiaL/measurementsTests"
hdf5path = ""
parameterfile = "input/parameters.csv"
main_write_hdf5_logs(path,hdf5path,parameterfile)