using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HDF5
using Plots
plotlyjs()

# In these correlators the periodicity is such that c[2] = c[T]
basepath = "./output/correlation_matrix/"
name = "Lt64Ls20beta6.5mf0.70mas1.01"
name = "Lt64Ls20beta6.5mf0.71mas1.01"
name = "Lt48Ls20beta6.5mf0.71mas1.01"

Γ = "g5"
file = joinpath(basepath,"correlation_matrix_$name.h5")
corr = h5read(file,"singlet_correlation_matrix_$Γ")

# crossing energy level at swap
ev, Δev, vecs, Δvecs = eigenvalues_eigenvectors(corr;t0=1,swap=nothing)
evJK, vecsJK = MixedRepSinglets.eigenvalues_eigenvectors_jackknife_samples(corr;t0=6)

# Pick a specific jackknife sample for testing
t = 8
sample = 1
vec1 = vecsJK[:,1,sample,t]
vec2 = vecsJK[:,2,sample,t]


vec1[1]/vec2[1]
vec2[2]/vec1[2]

# eigenvectors are still correctly normalized


# now we only need to fit the matrix at a given time to a rotation matrix
plot()
scatter!(vecs[1,1,:], yerr = Δvecs[1,1,:], label="(1,1)")
scatter!(vecs[2,1,:], yerr = Δvecs[2,1,:], label="(2,1)")
scatter!(vecs[1,2,:], yerr = Δvecs[1,2,:], label="(1,2)")
scatter!(vecs[2,2,:], yerr = Δvecs[2,2,:], label="(2,2)")
ylims!((-1.2,1.2))
xlims!((4,15))
