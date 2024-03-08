using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Plots
using Statistics
gr(legend=:bottomleft, frame=:box, legendfontsize=12, tickfontsize=12, labelfontsize=18, markersize=5)

Nsmear = 0:40:80
markershape = :rect
h5file_disc = "out_spectrum_smeared_discon_h128.hdf5"

Nsmear = 0:10:80
markershape = :circle
h5file_disc = "out_spectrum_smeared_discon_h128_old.hdf5"

h5data = "/home/fabian/Downloads/data.hdf5"
data_disc = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90/out_spectrum_discon/DISCON_SEMWALL SINGLET_g5_disc_re"

loop = [_get_disconnected_at_smearing_level(h5file_disc,N,"g5","FUN") for N in Nsmear]
disc = [4*unbiased_estimator(loop[i];subtract_vev=false,rescale=1) for i in eachindex(Nsmear)] 
disc = cat(disc...,dims=3)

loop_old = h5read(h5data,data_disc)
disc_old = unbiased_estimator(loop_old;subtract_vev=false,rescale=1)

d , Δd  = stdmean(disc,dims=1)
d0, Δd0 = stdmean(disc_old,dims=1)

#plt1 = plot(yscale=:log10,legend=:outerright)
scatter!(plt1,d0,yerr=Δd0,label="disc (0)";markershape)
for i in eachindex(Nsmear)
    scatter!(plt1,d[:,i] ,yerr=Δd[:,i] ,label="disc ($i)";markershape)
end
display(plt1)
