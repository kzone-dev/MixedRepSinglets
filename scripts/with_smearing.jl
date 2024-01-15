using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Plots
plotlyjs()

file = "/home/fabian/Documents/Lattice/HiRepDIaL/measurements/Lt64Ls20beta6.5mf0.71mas1.01FUN/out/out_spectrum_smeared"

N_max = 20
Nstep = 10
nhits = 8

h5fileCONN= "test_conn.hdf5"
h5fileDISC= "test_disc.hdf5"

function write_hdf5_smeared_combined(file,h5fileCONN,h5fileDISC,N_max,Nstep,nhits)
    typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in 0:Nstep:N_max]
    typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in 0:Nstep:N_max, N2 in 0:Nstep:N_max]
    writehdf5_spectrum(file,h5fileCONN,typesCONN)
    writehdf5_spectrum_disconnected(file,h5fileDISC,typesDISC,nhits)    
end
function _get_connected_at_smearing_level(h5file,Nsource,Nsink,channel)
    group = "source_N$(Nsource)_sink_N$(Nsink) TRIPLET"
    return h5read(h5file,joinpath(group,channel))
end
function _get_disconnected_at_smearing_level(h5file,Nsmear,channel)
    group = "DISCON_SEMWALL smear_N$Nsmear SINGLET"
    return h5read(h5file,joinpath(group,channel))
end
function _read_diagrams_smeared_single_repre(h5fileCONN,h5fileDISC,N1,N2;kws...)
    connN1N2 = _get_connected_at_smearing_level(h5fileCONN,N1,N2,channel)
    discN1 = _get_disconnected_at_smearing_level(h5fileDISC,N1,channel)
    discN2 = _get_disconnected_at_smearing_level(h5fileDISC,N2,channel)

    #rescale disconnected pieces to match the common normalisation
    T, L  = h5read(h5fileCONN,"lattice")[1:2]
    rescale_disc = (L^3)^2 /L^3
    if N1 == N2
        discN1N2 = unbiased_estimator(discN1;rescale=rescale_disc,kws...)
    else
        discN1N2 = unbiased_estimator(discN1,discN2;rescale=rescale_disc,kws...)
    end
    # rescale now the connected piece appropriately
    MixedRepSinglets.rescale_connected!(connN1N2,L)
    return connN1N2, discN1N2
end

N1 = 10
N2 = 20
channel = "g5"
conn, disc = _read_diagrams_smeared_single_repre(h5fileCONN,h5fileDISC,N1,N2)
