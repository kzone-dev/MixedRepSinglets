using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Plots
using Statistics
plotlyjs()

N_max  = 20
Nstep  = 10
Nsmear = 0:Nstep:N_max
nhits  = 32

fileCONN = "/home/fabian/Downloads/out_spectrum_smeared"
fileDISC = "/home/fabian/Downloads/out_spectrum_smeared_discon"
fileCONN_noAPE = "/home/fabian/Downloads/out_spectrum_smeared_noAPE"
fileDISC_noAPE = "/home/fabian/Downloads/out_spectrum_smeared_discon_noAPE"

h5file = "single_rep_smeared.hdf5"
h5file_noAPE = "single_rep_smeared_noAPE.hdf5"

typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]
typesDISC_noAPE = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in [0,10]]
typesCONN_noAPE = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in [0,10], N2 in [0,10]]

# add data from publication as reference
h5data = "/home/fabian/Downloads/data.hdf5"
group = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90"
file1 = "out_spectrum"
file2 = "out_spectrum_discon"
channel1 = "DEFAULT_SEMWALL TRIPLET_g5"
channel2 = "DISCON_SEMWALL SINGLET_g5_disc_re"

writehdf5 = false
if writehdf5
    writehdf5_spectrum(fileCONN,h5file,typesCONN,h5group="FUN/CONN",setup=true)
    writehdf5_spectrum_disconnected(fileDISC,h5file,typesDISC,nhits,h5group="FUN/DISC",setup=false)
    writehdf5_spectrum(fileCONN_noAPE,h5file_noAPE,typesCONN_noAPE,h5group="FUN/CONN",setup=true)
    writehdf5_spectrum_disconnected(fileDISC_noAPE,h5file_noAPE,typesDISC_noAPE,nhits,h5group="FUN/DISC",setup=false)
end

function _get_connected_at_smearing_level(h5file,Nsource,Nsink,channel,rep)
    group = "source_N$(Nsource)_sink_N$(Nsink) TRIPLET"
    return h5read(h5file,joinpath(rep,"CONN",group,channel))
end
function _get_disconnected_at_smearing_level(h5file,Nsmear,channel,rep)
    group = "DISCON_SEMWALL smear_N$Nsmear SINGLET"
    return h5read(h5file,joinpath(rep,"DISC",group,channel))
end
function _smeared_correlation_matrix(h5file,Nsmear,channel,rep;Nf_fun=2, disc_sign=+1, rescale_disc=1, rescale_conn = false)
    discFUN = [_get_disconnected_at_smearing_level(h5file,N,channel,rep) for N in Nsmear]
    connFUN = [_get_connected_at_smearing_level(h5file,N1,N2,channel,rep) for N1 in Nsmear, N2 in Nsmear ]

    # choose the smallest value of N for all measurements
    N1 = minimum(first.(size.(discFUN)))
    N2 = minimum(first.(size.(connFUN)))
    N  = minimum((N1,N2))
    T,L = h5read(h5file,"lattice")[1:2]

    # rescale connected
    rescale_conn && MixedRepSinglets.rescale_connected!.(connFUN,L)

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
    return corrMatCONN, corrMatDISC
end
function stdmean(X;dims)
    N = size(X)[dims]
    m = dropdims(mean(X;dims);dims)
    s = dropdims(std(X;dims);dims)/sqrt(N)
    return m, s
end

rescale_disc = 1
rescale_conn = false
T,L = h5read(h5file,"lattice")[1:2]

corrMatCONN, corrMatDISC = _smeared_correlation_matrix(h5file,Nsmear,"g5","FUN";rescale_disc,rescale_conn)
avgMatCONN, stdMatCONN = stdmean(corrMatCONN;dims=3)
avgMatDISC, stdMatDISC = stdmean(corrMatDISC;dims=3)

# add data from publication as reference
conn_old = h5read(h5data,joinpath(group,file1,channel1))
MixedRepSinglets.rescale_connected!(conn_old,1)
corr_conn, corr_conn_Delta = stdmean(conn_old;dims=2)

disc_old = h5read(h5data,joinpath(group,file2,channel2))
corr_disc_old = unbiased_estimator(disc_old;rescale=rescale_disc)
corr_disc, corr_disc_Delta = stdmean(corr_disc_old,dims=1)

# get data without APE smearing
discFUN_noAPE = _get_disconnected_at_smearing_level(h5file_noAPE,0,"g5","FUN")
connFUN_noAPE = _get_connected_at_smearing_level(h5file_noAPE,0,0,"g5","FUN") 
discFUN_corr_noAPE = unbiased_estimator(discFUN_noAPE;rescale=rescale_disc)
corr_disc_noAPE, corr_disc_noAPE_Delta = stdmean(discFUN_corr_noAPE,dims=1)
corr_conn_noAPE, corr_conn_noAPE_Delta = stdmean(connFUN_noAPE,dims=1)

plt1 = plot()
plt2 = plot()
for i in 1:3
    scatter!(plt1,0.5avgMatCONN[i,i,:],yerr=0.5stdMatCONN[i,i,:],label="N=$(Nsmear[i]) (with APE)(conn.)")
    scatter!(plt2,2avgMatDISC[i,i,:],yerr=2stdMatDISC[i,i,:],label="N=$(Nsmear[i]) (with APE)(disc.)")
end
plt1
scatter!(plt1,corr_conn,yerr=corr_conn_Delta, label="no APE no Wuppertal (old)")
scatter!(plt2,corr_disc,yerr=corr_disc_Delta, label="no APE no Wuppertal (old)")
scatter!(plt1,0.5corr_conn_noAPE,yerr=0.5corr_conn_noAPE_Delta, label="no APE no Wuppertal")
scatter!(plt2,4corr_disc_noAPE,yerr=4corr_disc_noAPE_Delta, label="no APE no Wuppertal")
plot!(plt1,yscale=:log10)
plot!(plt2,yscale=:log10)