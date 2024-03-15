using MixedRepSinglets
include("utils.jl")

hdf5path       = "/home/fabian/Downloads/hdf5out_modified"
parameters_gevp = joinpath(paramter_path,"parameters_gevp.csv")

h5corrs     = joinpath(hdf5path,"singlets_smeared_correlators.hdf5")
h5eigenvals = joinpath(hdf5path,"singlets_smeared_eigenvalues.hdf5")

parameters = readdlm(parameters_gevp,';';skipstart=1)
for row in eachrow(parameters)

    ensemble, channel, t0, binsize, deriv, ops = row
    nops = [1,10]
    deriv = false
    binsize = 2
    t0    = 1
    # get only entries without Wuppertal smearing

    channel == "g5_singlet" || continue
    
    matrixname ="correlation_matrix_g0g5_singlet"
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,matrixname))
    correlation_matrix = correlation_matrix[nops,nops,:,:]
    
    eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife, eigvecs, Δeigvecs = eigenvalues_eigenvectors_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    @show size(eigvals)

    β   = h5read(h5corrs,joinpath(ensemble,"beta"))
    T,L = h5read(h5corrs,joinpath(ensemble,"lattice"))[1:2]
    mf  = h5read(h5corrs,joinpath(ensemble,"quarkmasses_fundamental"))[1]
    mas = h5read(h5corrs,joinpath(ensemble,"quarkmasses_antisymmetric"))[1]

    plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=2,tag=L"$J^P = 0^-$(singlet)")
    plot!(plt1, ylims=(0.3,1.2))
    #display(plt1)

    asin_deriv(x) = +1/sqrt(1-x^2)
    acos_deriv(x) = -1/sqrt(1-x^2)
    asin_error(x,Δx) = abs(asin_deriv(x))*Δx
    acos_error(x,Δx) = abs(acos_deriv(x))*Δx

    @show size(eigvecs)
    plt3 = plot(legend=:outerright)
    scatter!(plt3, eigvecs[1,1,:],  ms = 8, yerr = Δeigvecs[1,1,:],label="(1,1)")
    scatter!(plt3, eigvecs[1,2,:],  ms = 8, yerr = Δeigvecs[1,2,:],label="(1,2)")
    scatter!(plt3, eigvecs[2,1,:],  ms = 5, yerr = Δeigvecs[2,1,:],label="(2,1)")
    scatter!(plt3, eigvecs[2,2,:],  ms = 5, yerr = Δeigvecs[2,2,:],label="(2,2)")
    plot!(xlims=(t0,10),ylims=(0.98,1.02))
    display(plt3)
end

   