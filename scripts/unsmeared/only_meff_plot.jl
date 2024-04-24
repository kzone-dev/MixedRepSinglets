using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HDF5
using Statistics
using Plots
using LaTeXStrings
#gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=11, tickfontsize=10, labelfontsize=14, markersize=6)

basepath  = "./output/correlation_matrix/"
ensembles = [
    "Lt48Ls20beta6.5mf0.71mas1.01",
    "Lt64Ls20beta6.5mf0.70mas1.01",
    "Lt64Ls20beta6.5mf0.71mas1.01",
    "Lt80Ls20beta6.5mf0.71mas1.01",
    "Lt96Ls20beta6.5mf0.71mas1.01"
]
swaps = [6, 6, 6, 6, 6]
tmax1 = [12, 12, 12, 12, 12] 
tmax2 = [21, 15, 15, 16, 18]

for (i,name) in enumerate(ensembles)
    file = joinpath(basepath,"correlation_matrix_$name.h5")
    corr = h5read(file,"singlet_correlation_matrix_g5")
    corr = h5read(file,"singlet_correlation_matrix_g5_folded_vev_sub")
    corr = h5read(file,"singlet_correlation_matrix_g5_folded")
    #corr = h5read(file,"singlet_correlation_matrix_id_folded_vev_sub")
    corr = _bin_correlator_matrix(corr;binsize=2) 
    sign = +1    
    #corr = correlator_derivative(corr;t_dim=4)
    #sign = -1

    swap = swaps[i]
    # obtain effective masses from jackknife analysis
    jks = eigenvalues_jackknife_samples(corr)
    m, Δm = meff_from_jackknife(jks;sign,swap)
    
    N = size(corr)[3]
    plt = plot(title=L" %$name ($N_{\rm eff}=%$N$)", xlabel="t", ylabel="effective mass")
    scatter!(plt,m[1,1:tmax2[i]], yerr=Δm[1,1:tmax2[i]],marker=:diamond,label=L"$m _{\rm eff}: a $")
    scatter!(plt,m[2,1:tmax1[i]], yerr=Δm[2,1:tmax1[i]],marker=:pentagon,label=L"$m _{\rm eff}: \eta'$")
    plot!(plt,ylims=(0.2,1))
    display(plt)

    # path for Plots
    path = "./plots/singlet_meff/"
    ispath(path) || mkpath(path)
    savefig(plt,joinpath(path,"$name.pdf"))
end