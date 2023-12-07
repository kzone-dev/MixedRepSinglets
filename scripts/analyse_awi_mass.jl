using MixedRepSinglets
using HDF5
using Plots
using Statistics
gr()

path = "./output/h5files/"
names = [ 
    "Lt48Ls20beta6.5mf0.71mas1.01AS",
    "Lt48Ls20beta6.5mf0.71mas1.01FUN",
    "Lt64Ls20beta6.5mf0.70mas1.01AS",
    "Lt64Ls20beta6.5mf0.70mas1.01FUN",
    "Lt64Ls20beta6.5mf0.71mas1.01AS",
    "Lt64Ls20beta6.5mf0.71mas1.01FUN",
    "Lt80Ls20beta6.5mf0.71mas1.01AS",
    "Lt80Ls20beta6.5mf0.71mas1.01FUN",
    "Lt96Ls20beta6.5mf0.71mas1.01AS",
    "Lt96Ls20beta6.5mf0.71mas1.01FUN"
]
tcut = [10, 8, 12, 8, 12, 10, 15, 10, 20, 20]

for (i,name) in enumerate(names)
    file = joinpath(path,name,"out_spectrum_mixed.h5")
    corr = awi_corr(file)
        
    N,T = size(corr)
    c  = dropdims(mean(corr,dims=1),dims=1)
    Δc = dropdims(std(corr,dims=1),dims=1)/sqrt(N)

    tmin= T÷2 - tcut[i]
    tmax= T÷2 + tcut[i] + 1

    m, Δm = MixedRepSinglets.awi_fit_jackknife(file;tmin,tmax,binsize=2)
    plt = scatter(c,yerr=Δc,label="data: AWI quark mass")
    plot!(plt,tmin:tmax,m*ones(tmax-tmin+1),ribbon=Δm,label="fit")
    plot!(plt,xlims=(5,T-5),ylims=(m - 5Δm,m + 5Δm))
    display(plt)
end