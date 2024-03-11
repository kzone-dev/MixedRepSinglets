using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
using LaTeXStrings
include("utils.jl")
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=11, tickfontsize=10, labelfontsize=14, markersize=5)
plotlyjs()

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M1234.hdf5"
h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M1234_with_conn.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M1234.hdf5"

function _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas)
    Nops = first(size(meff))

    rel_error0 = abs.(Δmeff[Nops-0,:]./meff[Nops-0,:])
    rel_error1 = abs.(Δmeff[Nops-1,:]./meff[Nops-1,:])
    tmax0 = findfirst(x-> x > (1/2), rel_error0)
    tmax1 = findfirst(x-> x > (1/2), rel_error1)
    range0 = 2:tmax0
    range1 = 2:tmax1

    title = L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"

    plt1 = plot(title=title, xlabel="t", ylabel="effective mass")
    scatter!(plt1,range0, meff[Nops-0,range0], yerr= Δmeff[Nops-0,range0],label=L"$m _{\rm eff}: a $")
    scatter!(plt1,range1, meff[Nops-1,range1], yerr= Δmeff[Nops-1,range1],label=L"$m _{\rm eff}: \eta'$")
    plot!(plt1, ylims=(0.2,1))

    plt2 = plot(title=title, xlabel="t", ylabel="eigenvalues")
    scatter!(plt2,eigvals[Nops-0,:], yerr= Δeigvals[Nops-0,:],label=L"$m _{\rm eff}: a $")
    scatter!(plt2,eigvals[Nops-1,:], yerr= Δeigvals[Nops-1,:],label=L"$m _{\rm eff}: \eta'$")

    return plt1, plt2
end

for ensemble in ["M1","M2","M3","M4"]
    
    t0      = 1
    binsize = 2
    deriv   = true

    write = false
    write && write_eigenvalues_and_effective_masses(h5eigenvals,h5corrs,ensemble,"correlation_matrix_g5_singlet";t0,binsize,deriv)
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(h5corrs,ensemble,"correlation_matrix_g5_singlet";t0,binsize,deriv)
    #eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(h5corrs,ensemble,"correlation_matrix_g5_nonsinglet_FUN";t0,binsize,deriv)
    #eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(h5corrs,ensemble,"correlation_matrix_g1_nonsinglet_FUN";t0,binsize,deriv)
    #eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(h5corrs,ensemble,"correlation_matrix_g5_nonsinglet_AS";t0,binsize,deriv)
    #eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(h5corrs,ensemble,"correlation_matrix_g1_nonsinglet_AS";t0,binsize,deriv)

    β   = h5read(h5corrs,joinpath(ensemble,"beta"))
    T,L = h5read(h5corrs,joinpath(ensemble,"lattice"))[1:2]
    mf  = h5read(h5corrs,joinpath(ensemble,"quarkmasses_fundamental"))[1]
    mas = h5read(h5corrs,joinpath(ensemble,"quarkmasses_antisymmetric"))[1]

    plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas)
    display(plt1)
    display(plt2)

    #isdir("plots/singlet_meff_smeared") || mkdir("plots/singlet_meff_smeared")
    #savefig(plt1,"plots/singlet_meff_smeared/$(ensemble).pdf")
end




