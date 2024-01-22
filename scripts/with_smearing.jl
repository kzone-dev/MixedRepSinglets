using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Plots
plotlyjs()

fileFUN = "/home/fabian/Documents/Lattice/HiRepDIaL/measurements/Lt64Ls20beta6.5mf0.71mas1.01AS/out/out_spectrum_smeared"
fileAS  = "/home/fabian/Documents/Lattice/HiRepDIaL/measurements/Lt64Ls20beta6.5mf0.71mas1.01FUN/out/out_spectrum_smeared"

N_max  = 20
Nstep  = 10
Nsmear = 0:Nstep:N_max
nhits  = 8

h5file = "smeared.hdf5"
writehdf = false

typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]

if writehdf
    writehdf5_spectrum(fileFUN,h5file,typesCONN,h5group="FUN/CONN",setup=true)
    writehdf5_spectrum(fileAS ,h5file,typesCONN,h5group="AS/CONN" ,setup=false)
    writehdf5_spectrum_disconnected(fileFUN,h5file,typesDISC,nhits,h5group="FUN/DISC",setup=false)
    writehdf5_spectrum_disconnected(fileAS ,h5file,typesDISC,nhits,h5group="AS/DISC" ,setup=false)    
end
function _get_connected_at_smearing_level(h5file,Nsource,Nsink,channel,rep)
    group = "source_N$(Nsource)_sink_N$(Nsink) TRIPLET"
    return h5read(h5file,joinpath(rep,"CONN",group,channel))
end
function _get_disconnected_at_smearing_level(h5file,Nsmear,channel,rep)
    group = "DISCON_SEMWALL smear_N$Nsmear SINGLET"
    return h5read(h5file,joinpath(rep,"DISC",group,channel))
end
#function correlation_matrix_smeared(h5file,Nsmear,channel)
#end
discFUN = [_get_disconnected_at_smearing_level(h5file,N,"g5","FUN") for N in Nsmear]
discAS  = [_get_disconnected_at_smearing_level(h5file,N,"g5","AS")  for N in Nsmear]
# first  index: source smearing 
# second index: sink   smearing
connFUN = [_get_connected_at_smearing_level(h5file,N1,N2,"g5","FUN") for N1 in Nsmear, N2 in Nsmear ]
connAS  = [_get_connected_at_smearing_level(h5file,N1,N2,"g5","AS") for N1 in Nsmear, N2 in Nsmear ]
# choose the smallest value of N for all measurements
N1 = minimum(first.(size.(discFUN)))
N2 = minimum(first.(size.(discAS)))
N3 = minimum(first.(size.(connFUN)))
N4 = minimum(first.(size.(connAS)))
N  = minimum((N1,N2,N3,N4))
T,L = h5read(h5file,"lattice")[1:2]

# Compared to the old code, there is another factor of 2 per loop missing
rescale_disc = 4*L^3
# rescale connected pieces
MixedRepSinglets.rescale_connected!.(connFUN ,L)
MixedRepSinglets.rescale_connected!.(connAS  ,L)
# create correlation matrix 
Nf_fun = 2
Nf_as  = 3
disc_sign = +1
subtract_vev = false
Nops = 2*length(Nsmear)

#TODO: Write everything in a single matrix
# create block matrices of the full correlation matrix
block_diag_FUN = zeros((Nops÷2,Nops÷2,N,T))
block_diag_AS  = zeros((Nops÷2,Nops÷2,N,T))
# Only the disconnected diagrams appear here. Since they  are symmetric under an interchange of the two individual
# diagrams, the block-off-diagonals are identical.
block_mixed    = zeros((Nops÷2,Nops÷2,N,T))

# assemble block matrices
for ind1 in eachindex(Nsmear)
    for ind2 in eachindex(Nsmear)
        if ind1 == ind2
            discFUN_N1N2 = unbiased_estimator(discFUN[ind1];rescale=rescale_disc,subtract_vev)
            discAS_N1N2  = unbiased_estimator(discAS[ind1] ;rescale=rescale_disc,subtract_vev)
        else
            discFUN_N1N2 = unbiased_estimator(discFUN[ind1],discFUN[ind2];rescale=rescale_disc,subtract_vev) 
            discAS_N1N2  = unbiased_estimator(discAS[ind1] ,discAS[ind2] ;rescale=rescale_disc,subtract_vev) 
        end
        discFUNAS_N1N2   = unbiased_estimator(discFUN[ind1],discAS[ind2] ;rescale=rescale_disc,subtract_vev) 
        block_diag_FUN[ind1,ind2,:,:] = connFUN[ind1,ind2][1:N,:] - Nf_fun*disc_sign*discFUN_N1N2[1:N,:]
        block_diag_AS[ind1,ind2,:,:]  = connAS[ind1,ind2][1:N,:]  - Nf_as *disc_sign*discAS_N1N2[1:N,:]
        block_mixed[ind1,ind2,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*discFUNAS_N1N2[1:N,:]
    end
end

block_row_1 = vcat(block_diag_FUN,block_mixed)
block_row_2 = vcat(block_mixed,block_diag_AS)
correlation_matrix = hcat(block_row_1,block_row_2)