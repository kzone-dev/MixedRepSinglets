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
function _smeared_singlet_correlation_matrix(h5file,Nsmear,channel,rep; rescale_disc=1, rescale_conn = false)
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
            corrMatDISC[ind1,ind2,:,:] = discFUN_N1N2[1:N,:]
        end
    end
    return corrMatCONN, corrMatDISC
end
function _smeared_singlet_correlator(h5file,N,channel,rep;kws...)
    # get data without APE smearing
    loop = _get_disconnected_at_smearing_level(h5file,N,channel,rep)
    conn = _get_connected_at_smearing_level(h5file,N,N,channel,rep) 
    disc = unbiased_estimator(loop;kws...)
    # rescale 
    MixedRepSinglets.rescale_connected!(conn,L)
    return conn, disc
end
function stdmean(X;dims)
    N = size(X)[dims]
    m = dropdims(mean(X;dims);dims)
    s = dropdims(std(X;dims);dims)/sqrt(N)
    return m, s
end

T,L = h5read(h5file,"lattice")[1:2]
rescale_disc = L^3
rescale_conn = true
Nf=2

conn_matrix, disc_matrix = _smeared_singlet_correlation_matrix(h5file,Nsmear,"g5","FUN";rescale_disc,rescale_conn)
conn_noAPE,  disc_noAPE  = _smeared_singlet_correlator(h5file_noAPE,0,"g5","FUN";rescale=rescale_disc)

# add data from publication as reference
h5data = "/home/fabian/Downloads/data.hdf5"
group  = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90"
file1  = "out_spectrum"
file2  = "out_spectrum_discon"
channel1 = "DEFAULT_SEMWALL TRIPLET_g5"
channel2 = "DISCON_SEMWALL SINGLET_g5_disc_re"

conn_old = h5read(h5data,joinpath(group,file1,channel1))
disc_old = h5read(h5data,joinpath(group,file2,channel2))
MixedRepSinglets.rescale_connected!(conn_old,L)
corr_disc_old = unbiased_estimator(disc_old;rescale=rescale_disc)
corr_conn, corr_conn_Delta = stdmean(conn_old;dims=2)
corr_disc, corr_disc_Delta = stdmean(corr_disc_old,dims=1)

# built full singlet correlator
# PART 1: Smeared correlators
#         Factor 2 is missing relative factor (not the Nf factor)
# PART 2: Smeared correlators without APE smearing (Nf included)
# PART 3: Old published data 
corr1 = conn_matrix   - 4Nf * disc_matrix
corr2 = conn_noAPE    - 4Nf * disc_noAPE
corr3 = conn_old'     -  Nf * corr_disc_old
# permute dimensions so that the first index corresponds to Euclidean time
#                            the second index refers to the Monte-Carlo time
corr1 = permutedims(corr1,(4,3,1,2))
corr2 = permutedims(corr2,(2,1))
corr3 = permutedims(corr3,(2,1))

c1, Δc1 = stdmean(corr1,dims=2)
c2, Δc2 = stdmean(corr2,dims=2) 
c3, Δc3 = stdmean(corr3,dims=2)

# Plot effective mass
meff1, Δmeff1 = implicit_meff_jackknife(corr1;sign=+1)
meff2, Δmeff2 = implicit_meff_jackknife(corr2;sign=+1)
meff3, Δmeff3 = implicit_meff_jackknife(corr3;sign=+1)

# Perform GEVP
# transform to the old inidces
corr1 = permutedims(corr1,(3,4,2,1))
jks = eigenvalues_jackknife_samples(corr1)
m, Δm = meff_from_jackknife(jks;sign=+1)

plt1 = plot()
for i in 1:3
    scatter!(plt1,c1[i,i,:],yerr=Δc1[i,i,:],label="N=$(Nsmear[i]) (with APE)(conn.)")
end
scatter!(plt1,c2,yerr=Δc2, label="no APE no Wuppertal (old)")
scatter!(plt1,c3,yerr=Δc3, label="no APE no Wuppertal")
plot!(plt1,yscale=:log10)

plt3 = plot()
scatter!(plt3,meff3 ,yerr=Δmeff3)
scatter!(plt3,m[3,:],yerr=Δm[3,:])
plot!(xlims=(0,12),ylims=(0.3,1.1))