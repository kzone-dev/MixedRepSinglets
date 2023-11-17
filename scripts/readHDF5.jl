using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HDF5
using Statistics
using Plots
using LaTeXStrings
gr(legend=:topright, frame=:box)
# In these correlators the periodicity is such that c[2] = c[T]
basepath = "./output/h5files/"

name = "Lt64Ls20beta6.5mf0.70mas1.01"
name = "Lt48Ls20beta6.5mf0.71mas1.01"
name = "Lt64Ls20beta6.5mf0.71mas1.01"

fileCONN_fun = joinpath(basepath,"$(name)FUN/out_spectrum.h5")
fileCONN_as  = joinpath(basepath,"$(name)AS/out_spectrum.h5")
fileDISC_fun = joinpath(basepath,"$(name)FUN/out_spectrum_discon.h5")
fileDISC_as  = joinpath(basepath,"$(name)AS/out_spectrum_discon.h5")

corr_deriv = correlator_derivative(corr;t_dim=4) 
ev, Δev = eigenvalues(corr)

# naive diagonal correlators
c1  = dropdims(mean(corr[1,1,:,:],dims=1),dims=1)
Δc1 = dropdims(std(corr[1,1,:,:],dims=1),dims=1)
c2  = dropdims(mean(corr[2,2,:,:],dims=1),dims=1)
Δc2 = dropdims(std(corr[2,2,:,:],dims=1),dims=1)

# obtain effective masses from jackknife analysis
jks = eigenvalues_jackknife_samples(corr)
m, Δm = meff_from_jackknife(jks;sign=-1)

plt1 = plot(yscale=:log10, ylabel=L"diagonal: $C(t)$" )
plot_correlator!(plt1,c1,Δc1;label="FUN")
plot_correlator!(plt1,c2,Δc2;label="AS")

plt2 = plot(yscale=:log10, ylabel=L"eigenvalues: $C(t)$" )
plot_correlator!(plt2,ev[1,:],Δev[1,:];label="EV1")
plot_correlator!(plt2,ev[2,:],Δev[2,:];label="EV2")

plt3 = plot(legendfontsize=10, ylabel="effective mass")
scatter!(plt3,m[1,:], yerr=Δm[1,:],label=L"$m _{\rm eff}: \eta^{a}$")
scatter!(plt3,m[2,:], yerr=Δm[2,:],label=L"$m _{\rm eff}: \eta^{b}$")

l = @layout [a; b; c]
s = (480, 3*200)
xlim=(1,16)

# choose limits for effective masses
skip = 3
ylim = extrema(filter(isfinite,m[1,skip:maximum(xlim)]))
plot!(plt3,ylims=(0.3,1.1), xlabel=L"t")

# meake full plot
plt = plot(plt1,plt2,plt3,layout=l,size=s,xlims=xlim,plot_title=name)
savefig(plt,name*".pdf")
display(plt)