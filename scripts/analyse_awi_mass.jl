using MixedRepSinglets
using HDF5
using Plots
plotlyjs()

path = "./output/h5files/"
name = "Lt48Ls20beta6.5mf0.71mas1.01AS"
file = joinpath(path,name,"out_spectrum_mixed.h5")

c, Δc = awi_corr(file)

tmin=17
tmax=35
mq, Δmq = awi_fit(c,Δc;tmin,tmax)
scatter(c,yerr=Δc)
plot!(tmin:tmax,mq*ones(tmax-tmin+1),ribbon=Δmq)
