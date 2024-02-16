using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Plots
using Statistics
#plotlyjs(frame=:box)
pgfplotsx(legend=:bottomleft, frame=:box, legendfontsize=12, tickfontsize=12, labelfontsize=18, markersize=5)
include("smearing_tools.jl")

N_max  = 80
Nstep  = 40
Nsmear = 0:Nstep:N_max
nhits  = 128

path = "/media/fabian/HDD#3/rsyncout/rsyncVSC/runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90/out"
path = "/home/fabian/Downloads/smearedVSC/"

fileCONN = joinpath(path,"out_spectrum_smeared")
fileDISC = joinpath(path,"out_spectrum_smeared_discon_h128")

h5file = "single_rep_smeared.hdf5"

typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]

writehdf5 = false
if writehdf5
    writehdf5_spectrum(fileCONN,h5file,typesCONN,h5group="FUN/CONN",setup=true)
    writehdf5_spectrum_disconnected(fileDISC,h5file,typesDISC,nhits,h5group="FUN/DISC",setup=false)
end

T,L = h5read(h5file,"lattice")[1:2]
rescale_disc = L^3
rescale_conn = true
subtract_vev = false
sign=+1
Nf=2
Ns=length(Nsmear)

# read data with and without APE smearing
conn_matrix, disc_matrix = _smeared_singlet_correlation_matrix(h5file,Nsmear,"g5","FUN";subtract_vev,rescale_disc,rescale_conn)
conn_matrix

# add data from publication as reference
h5data = "/home/fabian/Downloads/data.hdf5"
data_conn = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90/out_spectrum/DEFAULT_SEMWALL TRIPLET_g5"
data_disc = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90/out_spectrum_discon/DISCON_SEMWALL SINGLET_g5_disc_re"
conn_old = h5read(h5data,data_conn)
loop_old = h5read(h5data,data_disc)
disc_old = unbiased_estimator(loop_old;subtract_vev, rescale=rescale_disc)
MixedRepSinglets.rescale_connected!(conn_old,L)

# built full singlet correlator
corr1 = conn_matrix - 4Nf * disc_matrix
corr3 = conn_old'   -  Nf * disc_old

corr1 = correlator_folding(corr1;t_dim=4,sign)
corr3 = correlator_folding(corr3;t_dim=2,sign)
conn  = correlator_folding(conn_matrix;t_dim=2,sign)
disc  = 4Nf*correlator_folding(disc_matrix;t_dim=2,sign)

# symmetry sign of correlators
corr1 = correlator_derivative(corr1;t_dim=4)
corr3 = correlator_derivative(corr3;t_dim=2)
sign = -1

# permute dimensions so that the first index corresponds to Euclidean time
# and the second index refers to the Monte-Carlo time
corr1 = permutedims(corr1,(4,3,1,2))
disc  = permutedims(disc ,(4,3,1,2))
conn  = permutedims(conn ,(4,3,1,2))
corr3 = permutedims(corr3,(2,1))

c1, Δc1 = stdmean(corr1,dims=2)
c3, Δc3 = stdmean(corr3,dims=2)
c3C, Δc3C = stdmean(conn_old,dims=2)
c3D, Δc3D = stdmean(Nf*disc_old,dims=2)
cD, ΔcD = stdmean(disc,dims=2)
cC, ΔcC = stdmean(conn,dims=2)

plt1 = plot()
scatter!(plt1,cC[:,1,1],yerr=ΔcC[:,1,1],yscale=:log10)
scatter!(plt1,cC[:,2,2],yerr=ΔcC[:,2,2],yscale=:log10)
scatter!(plt1,cD[:,1,1],yerr=ΔcD[:,1,1],yscale=:log10)
scatter!(plt1,cD[:,2,2],yerr=ΔcD[:,2,2],yscale=:log10)
scatter!(plt1,c3       ,yerr=Δc3       ,yscale=:log10)
scatter!(plt1,c1[:,1,1],yerr=Δc1[:,1,1],yscale=:log10)
display(plt1)

# Plot effective mass
samples = eigenvalues_jackknife_samples(corr1,t0=1)
m, Δm = meff_from_jackknife(samples;sign)
meff1, Δmeff1 = implicit_meff_jackknife(corr1;sign)
meff3, Δmeff3 = implicit_meff_jackknife(corr3;sign)

xlabel = "Euclidean time"
ylabel = "effective mass"
title  = "pseudoscalar singlet: with numerical derivative"

plt2 = plot(;xlabel,ylabel,title)
plot!(plt2,xlims=(0,12),ylims=(0.4,0.9))
scatter!(plt2,m[Ns,1:10],yerr=Δm[Ns,1:10], label="GEVP N=0...$N_max (with APE)")
scatter!(meff1[:,1,1], yerr=Δmeff1[:,1,1], label="")
scatter!(meff1[:,2,2], yerr=Δmeff1[:,2,2], label="")
hspan!(plt2,[0.604,0.616],alpha=0.5,label="published result")
display(plt2)