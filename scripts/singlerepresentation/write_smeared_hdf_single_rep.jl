using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Statistics

path   = "/home/fabian/Downloads/"
h5file = "/home/fabian/Downloads/single_rep_smeared.hdf5"

Nsmear_conn = 0:40:80
Nsmear_disc = 0:40:80
nhits  = 128

name_conn = "out_spectrum_smeared"
name_disc = "out_spectrum_smeared_discon"
fileCONN = joinpath(path,name_conn)
fileDISC = joinpath(path,name_disc)

typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear_disc]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear_conn, N2 in Nsmear_conn]

writehdf5_spectrum(fileCONN,h5file,typesCONN,h5group="FUN/CONN")
writehdf5_spectrum_disconnected(fileDISC,h5file,typesDISC,nhits,h5group="FUN/DISC")
