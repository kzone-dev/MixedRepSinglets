using MixedRepSinglets
using HDF5
using Plots
plotlyjs()

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

for name in names
    file = joinpath(path,name,"out_spectrum_mixed.h5")
    c, Δc = awi_corr(file)

    T = length(c)

    tmin= 15
    tmax= T -  15

    mq, Δmq = awi_fit(c,Δc;tmin,tmax)
    plt = scatter(c,yerr=Δc)
    plot!(plt,tmin:tmax,mq*ones(tmax-tmin+1),ribbon=Δmq)
    display(plt)
end