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

io = Base.stdout
io = open("output/PCACmasses.csv","w")
println(io,"ensemble,T,L,beta,m0,mPCAC,ΔmPCAC")
for (i,name) in enumerate(names)
    file = joinpath(path,name,"out_spectrum_mixed.h5")

    Nf = contains(name,"FUN") ? 2 : 3
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
    #display(plt)


    T, L = h5read(file,"lattice")[1:2]
    beta = h5read(file,"beta")[]
    mass = h5read(file,"quarkmasses")[]

    println(io,"$name,$T,$L,$beta,$mass,$m,$Δm")
end
close(io)