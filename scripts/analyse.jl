using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
using LaTeXStrings
plotlyjs()
include("utils.jl")
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=11, tickfontsize=10, labelfontsize=14, markersize=5)

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M1.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M1.hdf5"

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M2.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M2.hdf5"

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M34.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M34.hdf5"

for ensemble in ["M4"]
    
    t0      = 1
    binsize = 2
    deriv   = true

    write = false
    write && write_eigenvalues_and_effective_masses(h5eigenvals,h5corrs,ensemble;t0binsize,deriv)
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(h5corrs,ensemble;t0,binsize,deriv)

    range1 = 2:8
    range0 = 2:11
    Nops   = 18

    β   = h5read(h5corrs,joinpath(ensemble,"beta"))
    T,L = h5read(h5corrs,joinpath(ensemble,"lattice"))[1:2]
    mf  = h5read(h5corrs,joinpath(ensemble,"quarkmasses"))[1]
    mas = "-1.01" #h5read(h5corrs,joinpath(ensemble,"quarkmasses"))[1]
    title = L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"

    plt = plot(title=title, xlabel="t", ylabel="effective mass")
    scatter!(plt,range0, meff[Nops-0,range0], yerr= Δmeff[Nops-0,range0],label=L"$m _{\rm eff}: a $")
    scatter!(plt,range1, meff[Nops-1,range1], yerr= Δmeff[Nops-1,range1],label=L"$m _{\rm eff}: \eta'$")
    plot!(plt, ylims=(0.2,1))
    display(plt)

    #isdir("plots/singlet_meff_smeared") || mkdir("plots/singlet_meff_smeared")
    #savefig(plt,"plots/singlet_meff_smeared/$(ensemble).pdf")

end




