using Pkg; Pkg.activate(".")
using HiRepParsing
using DelimitedFiles

path = "/home/fabian/Documents/Physics/Lattice/HiRepDIaL/measurements"
path = "/media/fabian/External SSD/DataDiaL/measurements"
path = "/home/fabian/Bilder/measurements"

input = readdlm("input/parameters_smeared.csv",';',skipstart=1)
h5file = "/home/fabian/Downloads/smeared_singlets_M34.hdf5"

for prm in eachrow(input)

    dir, patternDISC, patternCONN, nhits, Nmax, rep, name = prm

    fileCONN = joinpath(path,dir,"out/out_spectrum_smeared")
    fileDISC = joinpath(path,dir,"out/out_spectrum_smeared_discon")

    Nsmear = 0:10:Nmax
    typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
    typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]
    
    @show fileCONN   
    writehdf5_spectrum(fileCONN,h5file,typesCONN,h5group="$name/$rep/CONN")
    @show fileDISC   
    writehdf5_spectrum_disconnected(fileDISC,h5file,typesDISC,nhits,h5group="$name/$rep/DISC")
end    
