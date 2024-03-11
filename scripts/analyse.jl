using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
using LaTeXStrings
include("utils.jl")
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=11, tickfontsize=10, labelfontsize=14, markersize=5)

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M1234.hdf5"
h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M1234_with_conn.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M1234.hdf5"

function _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=1,tmax=nothing)
    Nops = first(size(meff))

    title = L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"   
    plt1 = plot(title=title, xlabel="t", ylabel="effective mass")
    plt2 = plot(title=title, xlabel="t", ylabel="eigenvalues")

    for state in 0:nstates-1
        if isnothing(tmax) 
            rel_error = abs.(Δmeff[Nops-state,:]./meff[Nops-state,:])
            tmax  = findfirst(x-> x > (1/2), rel_error)
        end
        range = 1:tmax
        scatter!(plt1,range, meff[Nops-state,range], yerr= Δmeff[Nops-state,range],label=L"$m _{\rm eff}$: state %$(state+1)")
        scatter!(plt2,range, eigvals[Nops-state,range], yerr= Δeigvals[Nops-state,range],label="EV: state $(state+1)")
    end
    plot!(plt1, ylims=(0.2,1))
    return plt1, plt2
end

for ensemble in ["M1"] #,"M2","M3","M4"]
    
    t0      = 1
    binsize = 2
    deriv   = true

    write = false
    write && write_eigenvalues_and_effective_masses(h5eigenvals,h5corrs,ensemble,"correlation_matrix_g5_singlet";t0,binsize,deriv)

    β   = h5read(h5corrs,joinpath(ensemble,"beta"))
    T,L = h5read(h5corrs,joinpath(ensemble,"lattice"))[1:2]
    mf  = h5read(h5corrs,joinpath(ensemble,"quarkmasses_fundamental"))[1]
    mas = h5read(h5corrs,joinpath(ensemble,"quarkmasses_antisymmetric"))[1]

    correlation_matrix = h5read(h5corrs,joinpath(ensemble,"correlation_matrix_g5_singlet"))
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=2)
    display(plt1)
    
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,"correlation_matrix_g5_nonsinglet_FUN"))
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv=false)
    plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=1,tmax=T÷2)
    display(plt1)
    
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,"correlation_matrix_g1_nonsinglet_FUN"))
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv=false)
    plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=1,tmax=T÷2)
    display(plt1)
    
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,"correlation_matrix_g5_nonsinglet_AS"))
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv=false)
    plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=1,tmax=T÷2)
    display(plt1)
    
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,"correlation_matrix_g1_nonsinglet_AS"))
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv=false)
    plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=1,tmax=T÷2)
    display(plt1)

    
    #isdir("plots/singlet_meff_smeared") || mkdir("plots/singlet_meff_smeared")
    #savefig(plt1,"plots/singlet_meff_smeared/$(ensemble).pdf")
end




