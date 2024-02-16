using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Plots
using Statistics
gr(legend=:bottomleft, frame=:box, legendfontsize=12, tickfontsize=12, labelfontsize=18, markersize=5)
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
h5data = "/home/fabian/Downloads/data.hdf5"
data_conn = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90/out_spectrum/DEFAULT_SEMWALL TRIPLET_g5"
data_disc = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90/out_spectrum_discon/DISCON_SEMWALL SINGLET_g5_disc_re"

typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]
writehdf5 = false
if writehdf5
    writehdf5_spectrum(fileCONN,h5file,typesCONN,h5group="FUN/CONN",setup=true)
    writehdf5_spectrum_disconnected(fileDISC,h5file,typesDISC,nhits,h5group="FUN/DISC",setup=false)
end

T,L = h5read(h5file,"lattice")[1:2]
sign=+1
Nf=2

# read data with and without APE smearing
conn, disc = _smeared_singlet_correlation_matrix(h5file,Nsmear,"g5","FUN";subtract_vev=false,rescale_disc=1,rescale_conn=false)
disc = 4Nf*disc

# add data from publication as reference
conn_old = h5read(h5data,data_conn)
loop_old = h5read(h5data,data_disc)
disc_old = unbiased_estimator(loop_old;subtract_vev=false,rescale=1)

# built full singlet correlator
corr1 = conn        - 4Nf * disc
corr3 = conn_old'   -  Nf * disc_old

c1,  Δc1  = stdmean(corr1,dims=3)
c3,  Δc3  = stdmean(corr3,dims=1)
c3C, Δc3C = stdmean(conn_old,dims=2)
c3D, Δc3D = stdmean(Nf*disc_old,dims=1)
cD,  ΔcD = stdmean(disc,dims=3)
cC,  ΔcC = stdmean(conn,dims=3)

plt1 = plot(legend=:outerright)
scatter!(plt1,cC[1,1,:],yerr=ΔcC[1,1,:],yscale=:log10,label="connected (1)")
scatter!(plt1,cC[2,2,:],yerr=ΔcC[2,2,:],yscale=:log10,label="connected (2)")
scatter!(plt1,cD[1,1,:],yerr=ΔcD[1,1,:],yscale=:log10,label="disconnected (1)")
scatter!(plt1,cD[2,2,:],yerr=ΔcD[2,2,:],yscale=:log10,label="disconnected (2)")
scatter!(plt1,c3C       ,yerr=Δc3C       ,yscale=:log10)
scatter!(plt1,c3D       ,yerr=Δc3D       ,yscale=:log10)
#scatter!(plt1,c3       ,yerr=Δc3       ,yscale=:log10)
#scatter!(plt1,c1[:,1,1],yerr=Δc1[:,1,1],yscale=:log10)
