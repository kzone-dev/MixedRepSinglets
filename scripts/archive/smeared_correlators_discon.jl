using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5
using Plots
using Statistics
#gr(legend=:bottomleft, frame=:box, legendfontsize=12, tickfontsize=12, labelfontsize=18, markersize=5)

Nsmear = 0:40:80
markershape = :rect
h5file_conn = "out_spectrum_smeared_more.hdf5"

Nsmear = 0:40:80
markershape = :circle
h5file_conn = "out_spectrum_smeared_single.hdf5"

Nsmear = 0:40:80
markershape = :pentagon
h5file_conn = "out_spectrum_smeared_git_51bc222fd.hdf5"

Nsmear = 0:40:80
markershape = :pentagon
h5file_conn = "out_spectrum_smeared_git_0bcb8163c.hdf5"

h5data = "/home/fabian/Downloads/data.hdf5"
data_conn = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90/out_spectrum/DEFAULT_SEMWALL TRIPLET_g5"

conn = [_get_connected_at_smearing_level(h5file_conn,N,N,"g5","FUN") for N in Nsmear]
conn = cat(conn...,dims=3)

conn_old = h5read(h5data,data_conn)

c , Δc  = stdmean(conn,dims=1)
c0, Δc0 = stdmean(conn_old,dims=2)

#plt1 = plot(yscale=:log10,legend=:outerright)
scatter!(plt1,c0,yerr=Δc0,label="conn (0)";markershape)
for i in eachindex(Nsmear)
    scatter!(plt1,c[:,i] ,yerr=Δc[:,i] ,label="conn ($i)";markershape)
end
display(plt1)
