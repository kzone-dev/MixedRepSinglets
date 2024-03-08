using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HDF5
using Statistics
using Plots
using LaTeXStrings
gr(legend=:topright, frame=:box, legendfontsize=12, tickfontsize=12, labelfontsize=14)

function overview_plot(c,Δc,ev,Δev,m,Δm;title="effective mass",s=(480, 3*200), xlim=(1,16),save=false, savedir="")
    c1,Δc1 = c[1,1,:], Δc[1,1,:]
    c2,Δc2 = c[2,2,:], Δc[2,2,:]
    
    plt1 = plot(yscale=:log10, ylabel=L"diagonal: $C(t)$" )
    plot_correlator!(plt1,c1,Δc1;label="FUN")
    plot_correlator!(plt1,c2,Δc2;label="AS")

    plt2 = plot(yscale=:log10, ylabel=L"eigenvalues: $C(t)$" )
    plot_correlator!(plt2,ev[1,:],Δev[1,:];label="EV1")
    plot_correlator!(plt2,ev[2,:],Δev[2,:];label="EV2")

    plt3 = plot(legendfontsize=10, ylabel="effective mass")
    scatter!(plt3,m[1,:], yerr=Δm[1,:],label=L"$m _{\rm eff}: \eta^{A}$")
    scatter!(plt3,m[2,:], yerr=Δm[2,:],label=L"$m _{\rm eff}: \eta^{B}$")

    l = @layout [a; b; c]

    # choose limits for effective masses
    skip = 3
    data = m[1,skip:maximum(xlim)]
    ylim = extrema(filter(x-> x >= 0 && isfinite(x),data))
    ylim = (0.25*ylim[1],1.33*ylim[2])
    plot!(plt3,ylims=ylim, xlabel=L"t")

    # meake full plot
    plt = plot(plt1,plt2,plt3,layout=l,size=s,xlims=xlim,plot_title=title)
    if save
        ispath(savedir) || mkpath(savedir)
        savefig(plt,joinpath(savedir,title*".pdf"))
    end
    display(plt)
    return plt
end

# In these correlators the periodicity is such that c[2] = c[T]
basepath = "./output/correlation_matrix/"

ensembles = [
    "Lt64Ls20beta6.5mf0.71mas1.01",
    "Lt64Ls20beta6.5mf0.70mas1.01",
    "Lt48Ls20beta6.5mf0.71mas1.01",
    "Lt80Ls20beta6.5mf0.71mas1.01",
    "Lt96Ls20beta6.5mf0.71mas1.01"
]
swaps = [nothing, nothing, nothing, nothing, nothing]

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

    overview_plot(c,Δc,ev,Δev,m,Δm,title=name,save=false,savedir="plots/correlators_effectivmasses/")
end