function all_effective_mass_plots(hdf5path,gevp_parameterfile)
    h5eigenvals = joinpath(hdf5path,"singlets_smeared_eigenvalues.hdf5")

    parameters = readdlm(gevp_parameterfile,';';skipstart=1)
    nrows = first(size(parameters))
    for row in 1:nrows

        #row > 4 && continue

        ensemble, channel, t0, binsize, deriv, ops = parameters[row,:]
        
        nops = parse.(Int,split(replace(ops,r"[()]"=>""),','))

        β   = h5read(h5eigenvals,joinpath(ensemble,channel,"beta"))
        T,L = h5read(h5eigenvals,joinpath(ensemble,channel,"lattice"))[1:2]
        mf  = h5read(h5eigenvals,joinpath(ensemble,channel,"quarkmasses_fundamental"))[1]
        mas = h5read(h5eigenvals,joinpath(ensemble,channel,"quarkmasses_antisymmetric"))[1]
        meff     = h5read(h5eigenvals,joinpath(ensemble,channel,"meff"))
        eigvals  = h5read(h5eigenvals,joinpath(ensemble,channel,"eigvals"))
        Δmeff    = h5read(h5eigenvals,joinpath(ensemble,channel,"Delta_meff"))
        Δeigvals = h5read(h5eigenvals,joinpath(ensemble,channel,"Delta_eigvals"))

        title = L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"   
        title = "ensemble $ensemble"
        plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,channel;title,nstates=2)
        plot!(plt1, ylims=(0.3,1.2))
        plot!(plt2,yscale=:log10)
        display(plt1)
        #display(plt2)
    end
end