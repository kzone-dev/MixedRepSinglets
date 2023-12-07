using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HDF5
using Statistics
using Plots
using LaTeXStrings
pgfplotsx(legend=:topright, frame=:box, legendfontsize=14, tickfontsize=12, labelfontsize=16, markersize=6)


# In these correlators the periodicity is such that c[2] = c[T]
basepath = "./output/correlation_matrix/"

ensembles = [
    "Lt64Ls20beta6.5mf0.71mas1.01",
    "Lt64Ls20beta6.5mf0.70mas1.01",
    "Lt48Ls20beta6.5mf0.71mas1.01",
    "Lt80Ls20beta6.5mf0.71mas1.01",
    "Lt96Ls20beta6.5mf0.71mas1.01"
]
swaps = [6, 6, 6, 6, 6]
tmax1 = [12, 12, 12, 12, 12] 
tmax2 = [21, 15, 15, 16, 18]

for (i,name) in enumerate(ensembles)
    file = joinpath(basepath,"correlation_matrix_$name.h5")
    corr = h5read(file,"singlet_correlation_matrix_g5")
    corr = h5read(file,"singlet_correlation_matrix_g5_folded")
    corr = _bin_correlator_matrix(corr;binsize=2) 

    swap = swaps[i]
    ev, Δev = eigenvalues(corr;swap)

    # naive diagonal correlators
    c  = dropdims(mean(corr,dims=3),dims=3)
    Δc = dropdims(std(corr,dims=3),dims=3)

    # obtain effective masses from jackknife analysis
    jks = eigenvalues_jackknife_samples(corr)
    m, Δm = meff_from_jackknife(jks;sign=+1,swap)

    plt3 = plot(title=name, xlabel="t", ylabel="effective mass")
    scatter!(plt3,m[1,1:tmax2[i]], yerr=Δm[1,1:tmax2[i]],marker=:diamond,label=L"$m _{\rm eff}: a $")
    scatter!(plt3,m[2,1:tmax1[i]], yerr=Δm[2,1:tmax1[i]],marker=:pentagon,label=L"$m _{\rm eff}: \eta'$")
    plot!(plt3,ylims=(0.2,1))

    display(plt3)
end