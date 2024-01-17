using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Plots
plotlyjs()

fileCONN = "/home/fabian/Downloads/out_spectrum_smeared"
fileDISC = "/home/fabian/Downloads/out_spectrum_smeared_discon"

N_max  = 20
Nstep  = 10
Nsmear = 0:Nstep:N_max
nhits  = 32

h5file = "single_rep_smeared.hdf5"

typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]
#writehdf5_spectrum(fileCONN,h5file,typesCONN,h5group="FUN/CONN",setup=true)
#writehdf5_spectrum_disconnected(fileDISC,h5file,typesDISC,nhits,h5group="FUN/DISC",setup=false)

function _get_connected_at_smearing_level(h5file,Nsource,Nsink,channel,rep)
    group = "source_N$(Nsource)_sink_N$(Nsink) TRIPLET"
    return h5read(h5file,joinpath(rep,"CONN",group,channel))
end
function _get_disconnected_at_smearing_level(h5file,Nsmear,channel,rep)
    group = "DISCON_SEMWALL smear_N$Nsmear SINGLET"
    return h5read(h5file,joinpath(rep,"DISC",group,channel))
end
discFUN = [_get_disconnected_at_smearing_level(h5file,N,"g5","FUN") for N in Nsmear]
connFUN = [_get_connected_at_smearing_level(h5file,N1,N2,"g5","FUN") for N1 in Nsmear, N2 in Nsmear ]

# choose the smallest value of N for all measurements
N1 = minimum(first.(size.(discFUN)))
N2 = minimum(first.(size.(connFUN)))
N  = minimum((N1,N2))
T,L = h5read(h5file,"lattice")[1:2]

# rescaling 
rescale_disc = 1  
MixedRepSinglets.rescale_connected!.(connFUN ,1)

# create correlation matrix 
Nf_fun = 2
disc_sign = +1
Nops = length(Nsmear)
# create block matrices of the full correlation matrix
corrMatCONN = zeros((Nops,Nops,N,T))
corrMatDISC = zeros((Nops,Nops,N,T))

# assemble block matrices
for ind1 in eachindex(Nsmear)
    for ind2 in eachindex(Nsmear)
        if ind1 == ind2
            discFUN_N1N2 = unbiased_estimator(discFUN[ind1];rescale=rescale_disc)
        else
            discFUN_N1N2 = unbiased_estimator(discFUN[ind1],discFUN[ind2];rescale=rescale_disc) 
        end
        corrMatCONN[ind1,ind2,:,:] = connFUN[ind1,ind2][1:N,:]
        corrMatDISC[ind1,ind2,:,:] = Nf_fun*disc_sign*discFUN_N1N2[1:N,:]
    end
end

using Statistics
avgMatCONN = dropdims(mean(corrMatCONN,dims=3),dims=3)
avgMatDISC = dropdims(mean(corrMatDISC,dims=3),dims=3)
stdMatCONN = dropdims(std(corrMatCONN,dims=3),dims=3)/sqrt(N)
stdMatDISC = dropdims(std(corrMatDISC,dims=3),dims=3)/sqrt(N)

using Plots
plt = plot()
for i in 1:3
    scatter!(plt,avgMatCONN[i,i,:],yerr=stdMatCONN[i,i,:],label="N=$(Nsmear[i])(conn.)")
    scatter!(plt,avgMatDISC[i,i,:],yerr=stdMatDISC[i,i,:],label="N=$(Nsmear[i])(disc.)")
end
plot!(plt,yscale=:log10)

# add unsmeared data for reference
h5data = "/home/fabian/Downloads/data.hdf5"
group = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90"
file1 = "out_spectrum"
file2 = "out_spectrum_discon"
channel1 = "DEFAULT_SEMWALL TRIPLET_g5"
channel2 = "DISCON_SEMWALL SINGLET_g5_disc_re"

conn_old = h5read(h5data,joinpath(group,file1,channel1))
N0 = size(conn_old)[2]
corr_conn = dropdims(mean(conn_old,dims=2),dims=2)
corr_conn_Delta = dropdims(std(conn_old,dims=2),dims=2)/sqrt(N0)

disc_old = h5read(h5data,joinpath(group,file2,channel2))
rescale_disc_old = 1
corr_disc_old = unbiased_estimator(disc_old;rescale=rescale_disc_old)
corr_disc = dropdims(mean(corr_disc_old,dims=1),dims=1)
corr_disc_Delta = dropdims(std(corr_disc_old,dims=1),dims=1)/sqrt(N0)

scatter!(plt,corr_conn,yerr=corr_conn_Delta, label="no APE")
scatter!(plt,corr_disc,yerr=2corr_disc_Delta, label="no APE")
