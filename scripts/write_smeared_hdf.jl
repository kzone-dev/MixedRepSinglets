using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Statistics
path = "/home/fabian/Downloads/smearedVSC/"
include("smearing_tools.jl")

Nsmear_conn = 0:10:10
Nsmear_disc = 0:10:80
nhits  = 128

name_conn = "out_spectrum_smeared_single_N10"
name_disc = "out_spectrum_smeared_discon_h128_old"
fileCONN = joinpath(path,name_conn)
fileDISC = joinpath(path,name_disc)

h5file_conn = "$name_conn.hdf5"
h5file_disc = "$name_disc.hdf5"

typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear_disc]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear_conn, N2 in Nsmear_conn]

#writehdf5_spectrum(fileCONN,h5file_conn,typesCONN,h5group="FUN/CONN")
writehdf5_spectrum_disconnected(fileDISC,h5file_disc,typesDISC,nhits,h5group="FUN/DISC")
