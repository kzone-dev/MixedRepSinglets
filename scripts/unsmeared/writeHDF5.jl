using Pkg; Pkg.activate(".")
using HiRepParsing
using DelimitedFiles

input = readdlm("./input/parameters.csv",';',skipstart=1)
filenames, typeDISC, typeCONN, nhits = input[:,1],input[:,2],input[:,3],input[:,4]

for (i,filename) in enumerate(filenames)
    basepath = "/home/fabian/Documents/Lattice/HiRepDIaL/measurements/"
    path   = joinpath(basepath,filename,"out")
    h5path = joinpath("./output/h5files",filename)

    ispath(h5path) || mkpath(h5path)

    # input files
    file_spectrum1 = joinpath(path,"out_spectrum_mixed")
    file_spectrum2 = joinpath(path,"out_spectrum")
    file_discon1  = joinpath(path,"out_spectrum_discon")
    
    # output files
    h5file_spectrum1 = joinpath(h5path,"out_spectrum_mixed.h5")
    h5file_spectrum2 = joinpath(h5path,"out_spectrum.h5")
    h5file_discon1  = joinpath(h5path,"out_spectrum_discon.h5")
    
    @show file_spectrum1
    writehdf5_spectrum(file_spectrum1,h5file_spectrum1,typeCONN[i])
    @show file_spectrum2
    writehdf5_spectrum(file_spectrum2,h5file_spectrum2,typeCONN[i])
    @show file_discon1
    writehdf5_spectrum_disconnected(file_discon1,h5file_discon1,typeDISC[i],nhits[i])
end