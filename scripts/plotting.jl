function plot_all_masses_with_fitting(parameters_gevp,parameters_fitting,corrfitterpath,hdf5path,plotdir;only_singlet=true)

    h5eigenvals = joinpath(hdf5path,"singlets_smeared_eigenvalues.hdf5")
    results_corrfitter = joinpath(corrfitterpath,"corrfitter_results.csv")

    dir1 = joinpath(plotdir,"effective_mass")
    dir2 = joinpath(plotdir,"effective_mass_groundstate")
    dir3 = joinpath(plotdir,"correlator")
    ispath(dir1) || mkpath(dir1)
    ispath(dir2) || mkpath(dir2)
    ispath(dir3) || mkpath(dir3)

    parameters = readdlm(parameters_gevp,';';skipstart=1)
    parameters_fitting = readdlm(parameters_fitting,';';skipstart=1)
    corrfitter_results = readdlm(results_corrfitter,';';skipstart=0)
    
    #check that the number of datasets match
    @assert first(size(parameters)) == first(size(parameters_fitting)) == first(size(corrfitter_results)) 
    nrows = first(size(parameters))

    function channel_tags(channel)
        isequal(channel,"g5_singlet")        && return L"$J^P = 0^-$(singlet)"
        isequal(channel,"g5_nonsinglet_FUN") && return L"$J^P = 0^-$(f)"
        isequal(channel,"g5_nonsinglet_AS")  && return L"$J^P = 0^-$(as)"
        isequal(channel,"g1_nonsinglet_FUN") && return L"$J^P = 1^+$(f)"
        isequal(channel,"g1_nonsinglet_AS")  && return L"$J^P = 1^+$(as)"
        return ""
    end

    for row in 1:nrows

        ensemble, channel, t0, binsize, deriv, ops = parameters[row,:]
        ensemble, channel, tmin1, tmin2, tmax1, tmax2, tp, Nmax  = parameters_fitting[row,:]
        ensemble, channel, T, L, mf, mas, beta, E0, ΔE0, E1, ΔE1, χ2dof0, χ2dof0  = corrfitter_results[row,:]

        nops = parse.(Int,split(replace(ops,r"[()]"=>""),','))

        β   = h5read(h5eigenvals,joinpath(ensemble,channel,"beta"))
        T,L = h5read(h5eigenvals,joinpath(ensemble,channel,"lattice"))[1:2]
        mf  = h5read(h5eigenvals,joinpath(ensemble,channel,"quarkmasses_fundamental"))[1]
        mas = h5read(h5eigenvals,joinpath(ensemble,channel,"quarkmasses_antisymmetric"))[1]
        meff     = h5read(h5eigenvals,joinpath(ensemble,channel,"meff"))
        eigvals  = h5read(h5eigenvals,joinpath(ensemble,channel,"eigvals"))
        Δmeff    = h5read(h5eigenvals,joinpath(ensemble,channel,"Delta_meff"))
        Δeigvals = h5read(h5eigenvals,joinpath(ensemble,channel,"Delta_eigvals"))

        title = "ensemble $ensemble" #L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"   
        plt1, plt2 = _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,channel;title,nstates=2)
        add_fit_range!(plt1,tmin1,tmax1,E0,ΔE0;label="")
        add_fit_range!(plt1,tmin2,tmax2,E1,ΔE1;label="")

        plot!(plt1, ylims=(0.9*E0,1.1*E0))
        #savefig(plt1,joinpath(dir2,"$(ensemble)_$(channel).pdf"))
        plot!(plt1, ylims=(0.3,1.2))
        isequal(channel,"g5_singlet") && savefig(plt1,joinpath(dir1,"$(ensemble)_$(channel).pdf"))
        plot!(plt2,yscale=:log10)
        #savefig(plt2,joinpath(dir3,"$(ensemble)_$(channel).pdf"))
        if !only_singlet || channel == "g5_singlet"
            display(plt1)
        end
    end
end
