function all_effective_mass_plots(h5eigenvals,gevp_parameterfile)

    parameters = readdlm(gevp_parameterfile,';';skipstart=1)
    nrows = first(size(parameters))
    for row in 1:nrows

        row > 4 && continue

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

        plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=2,tag=L"$J^P = 0^-$(singlet)")
        plot!(plt1, ylims=(0.3,1.2))
        plot!(plt2,yscale=:log10)
        display(plt1)
        display(plt2)
    end
end

