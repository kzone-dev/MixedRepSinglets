using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HDF5
using Plots
plotlyjs()

basepath = "./output/correlation_matrix/"
name = "Lt64Ls20beta6.5mf0.70mas1.01"
name = "Lt48Ls20beta6.5mf0.71mas1.01"
name = "Lt64Ls20beta6.5mf0.71mas1.01"

Γ = "g5"
file = joinpath(basepath,"correlation_matrix_$name.h5")
corr = h5read(file,"singlet_correlation_matrix_$Γ")

# crossing energy level at swap
t0   = 1
swap = nothing
evJK, vecsJK = MixedRepSinglets.eigenvalues_eigenvectors_jackknife_samples(corr;t0)

nops, nconf, T = size(evJK)
α1JK = zeros(nconf,T)
α2JK = zeros(nconf,T)
α3JK = zeros(nconf,T)

for t in 1:T, sample in 1:nconf
    vec1 = vecsJK[:,1,sample,t]
    vec2 = vecsJK[:,2,sample,t]
    # Three different definitions of the mixing angle
    # See e.g. eq. (17)-(20) in 1112.4384
    ratio1 = vec2[2]/vec1[2]
    ratio2 = -vec1[1]/vec2[1]

    α1JK[sample,t] = atand(ratio1)
    α2JK[sample,t] = atand(ratio2)
    α3JK[sample,t] = atand(sqrt(abs(ratio1*ratio2)))
end

α1, Δα1 = MixedRepSinglets.apply_jackknife(α1JK,dims=1)
α2, Δα2 = MixedRepSinglets.apply_jackknife(α2JK,dims=1)
α3, Δα3 = MixedRepSinglets.apply_jackknife(α3JK,dims=1)

plot()
scatter!(abs.(α1),yerr=Δα1,label="angle FUN")
scatter!(abs.(α2),yerr=Δα2,label="angle AS")
scatter!(abs.(α3),yerr=Δα3,label="angle AVG")
plot!(xlims=(6.5,13.5),ylims=(0,30))

#=
ev, Δev, vecs, Δvecs = eigenvalues_eigenvectors(corr;t0,swap)
plot()
scatter!(vecs[1,1,:], yerr = Δvecs[1,1,:], label="(1,1)")
scatter!(vecs[2,1,:], yerr = Δvecs[2,1,:], label="(2,1)")
scatter!(vecs[1,2,:], yerr = Δvecs[1,2,:], label="(1,2)")
scatter!(vecs[2,2,:], yerr = Δvecs[2,2,:], label="(2,2)")
ylims!((-1.2,1.2))
xlims!((4,15))
=#