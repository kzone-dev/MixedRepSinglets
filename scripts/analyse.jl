using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
using LaTeXStrings
plotlyjs()
include("utils.jl")
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=11, tickfontsize=10, labelfontsize=14, markersize=5)

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues.hdf5"

for ensemble in ["M1","M2","M3","M4","M5"]
    
    t0      = 1
    binsize = 2
    deriv   = true
    Nops   = 18

    write = true
    write && write_eigenvalues_and_effective_masses(h5eigenvals,h5corrs,ensemble;t0,binsize,deriv)
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(h5corrs,ensemble;t0,binsize,deriv)

    rel_error0 = abs.(Δmeff[Nops-0,:]./meff[Nops-0,:])
    rel_error1 = abs.(Δmeff[Nops-1,:]./meff[Nops-1,:])
    tmax0 = findfirst(x-> x > (1/2), rel_error0)
    tmax1 = findfirst(x-> x > (1/2), rel_error1)
    range0 = 2:tmax0
    range1 = 2:tmax1
 
    β   = h5read(h5corrs,joinpath(ensemble,"beta"))
    T,L = h5read(h5corrs,joinpath(ensemble,"lattice"))[1:2]
    mf  = h5read(h5corrs,joinpath(ensemble,"quarkmasses_fundamental"))[1]
    mas = h5read(h5corrs,joinpath(ensemble,"quarkmasses_antisymmetric"))[1]
    title = L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"

    plt = plot(title=title, xlabel="t", ylabel="effective mass")
    scatter!(plt,range0, meff[Nops-0,range0], yerr= Δmeff[Nops-0,range0],label=L"$m _{\rm eff}: a $")
    scatter!(plt,range1, meff[Nops-1,range1], yerr= Δmeff[Nops-1,range1],label=L"$m _{\rm eff}: \eta'$")
    plot!(plt, ylims=(0.2,1))
    display(plt)

    isdir("plots/singlet_meff_smeared") || mkdir("plots/singlet_meff_smeared")
    savefig(plt,"plots/singlet_meff_smeared/$(ensemble).pdf")

end




