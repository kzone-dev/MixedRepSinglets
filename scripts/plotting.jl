function plot_all_masses_with_fitting(parameters_gevp,parameters_fitting,corrfitterpath,hdf5path,plotdir;only_singlet=true)

    h5eigenvals = joinpath(hdf5path,"singlets_smeared_eigenvalues.hdf5")
    results_corrfitter = joinpath(corrfitterpath,"corrfitter_results.csv")

    dir1 = joinpath(plotdir,"effective_mass")
    dir2 = joinpath(plotdir,"eigvals")
    ispath(dir1) || mkpath(dir1)
    ispath(dir2) || mkpath(dir2)

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
        plot!(plt2,yscale=:log10)

        if !only_singlet || channel == "g5_singlet"
            savefig(plt1,joinpath(dir1,"$(ensemble)_$(channel).pdf"))
            savefig(plt2,joinpath(dir2,"$(ensemble)_$(channel).pdf"))
        end
    end
end
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
    end
end
function plot_spectrum(tablepath,plotdir,gradient_flow_results)
    results_MR = joinpath(tablepath,"table_results_MR.csv")

    dir = joinpath(plotdir,"spectrum")
    ispath(dir) || mkpath(dir)

    colours = palette(:Paired_6)

    # get meson masses
    data = readdlm(results_MR,';',skipstart=1)
    gf = readdlm(gradient_flow_results,',',skipstart=1)

    plt0 = plot(legend_column=2,legendfontsize=16,legend=:left,ylabel=L"am", xlabel=L"$am_0^{\rm f}$", title="meson masses vs. bare fermion mass")
    plt1 = plot(legend_column=2,legendfontsize=16,legend=:left,ylabel=L"am", xlabel=L"N_t", title=L"meson masses vs. spatial extent $N_t$")
    plt2 = plot(legend_column=2,legendfontsize=16,legend=:left,ylabel=L"am", xlabel=L"m_{\rm PS} / m_{\rm V}", title = L"meson masses vs. $m_{\rm PS} / m_{\rm V} ~ (N_t=64)$" )
    plt3 = plot(legend_column=2,legendfontsize=16,legend=:left,ylabel=L"wm", xlabel=L"m_{\rm PS} / m_{\rm V}", title=L"meson mass in gradient flow scale $w ~ (N_t=64)$")
    
    xticks0 = Float64[]
    xticks1 = Int[]

    for row in eachrow(data)
        ensemble, channel, T, L, mf, mas, ma, Δma, mη, Δmη, mπF, ΔmπF, mπA, ΔmπA, mρF, ΔmρF, mρA, ΔmρA = row

        ensemble_id = findfirst(isequal(ensemble), gf[:,1])
        w, Δw = gf[ensemble_id,7], gf[ensemble_id,8]

        Δratio(x,y,Δx,Δy) = sqrt((Δx/y)^2 + (x*Δy/y^2)^2) 
        Δproduct(x,y,Δx,Δy) = sqrt((Δx*y)^2 + (x*Δy)^2)
        r  = mπF/mρF
        Δr = Δratio(mπF,mρF,ΔmπF,ΔmρF)

        offset = sign(T-64)*0.0005

        push!(xticks0,mf)
        unique!(xticks0)

        push!(xticks1,T)
        unique!(xticks1)
 
        # add label with M2 ensemble, since it features in every plot
        label1 = (T == 64) && isapprox(mf,-0.71) ? L"\eta^{\prime}_h~~" : ""
        label2 = (T == 64) && isapprox(mf,-0.71) ? L"\eta^{\prime}_l~~" : ""
        label3 = (T == 64) && isapprox(mf,-0.71) ? L"{\rm ps}~~" : ""
        label4 = (T == 64) && isapprox(mf,-0.71) ? L"{\rm PS}~~" : ""
        label5 = (T == 64) && isapprox(mf,-0.71) ? L"{\rm v}~~" : ""
        label6 = (T == 64) && isapprox(mf,-0.71) ? L"{\rm V}~~" : ""
        
        scatter!(plt0, [mf+offset], [mη], yerr=Δmη, markersize=5,marker=:rect, color=colours[1];label=label1)
        scatter!(plt0, [mf+offset], [ma], yerr=Δma, markersize=5,marker=:rect, color=colours[2];label=label2)
        scatter!(plt0, [mf+offset], [mπA], yerr=ΔmπA, markersize=5,marker=:circ, color=colours[3];label=label3)
        scatter!(plt0, [mf+offset], [mπF], yerr=ΔmπF, markersize=5,marker=:circ, color=colours[4];label=label4)
        scatter!(plt0, [mf+offset], [mρA], yerr=ΔmρA, markersize=5,marker=:pentagon, color=colours[5];label=label5)
        scatter!(plt0, [mf+offset], [mρF], yerr=ΔmρF, markersize=5,marker=:pentagon, color=colours[6];label=label6)
        plot!(plt0;xticks=xticks0)

        if isapprox(mf,-0.71)
            scatter!(plt1, [T], [mη], yerr=Δmη, color=colours[1], markersize=5, marker=:rect;label=label1)
            scatter!(plt1, [T], [ma], yerr=Δma, color=colours[2], markersize=5, marker=:rect;label=label2)
            scatter!(plt1, [T], [mπA], yerr=ΔmπA, color=colours[3], markersize=5, marker=:circ;label=label3)
            scatter!(plt1, [T], [mπF], yerr=ΔmπF, color=colours[4], markersize=5, marker=:circ;label=label4)
            scatter!(plt1, [T], [mρA], yerr=ΔmρA, color=colours[5], markersize=5, marker=:pentagon;label=label5)
            scatter!(plt1, [T], [mρF], yerr=ΔmρF, color=colours[6], markersize=5, marker=:pentagon;label=label6)
            plot!(plt1;xticks=xticks1)
        end

        if T == 64
            scatter!(plt2, [r], xerr = Δr, [mη], yerr=Δmη, color=colours[1] ,markersize=5, marker=:rect;label=label1)
            scatter!(plt2, [r], xerr = Δr, [ma], yerr=Δma, color=colours[2] ,markersize=5, marker=:rect;label=label2)
            scatter!(plt2, [r], xerr = Δr, [mπA], yerr=ΔmπA, color=colours[3], markersize=5, marker=:circ;label=label3)
            scatter!(plt2, [r], xerr = Δr, [mπF], yerr=ΔmπF, color=colours[4], markersize=5, marker=:circ;label=label4)
            scatter!(plt2, [r], xerr = Δr, [mρA], yerr=ΔmρA, color=colours[5], markersize=5, marker=:pentagon;label=label5)
            scatter!(plt2, [r], xerr = Δr, [mρF], yerr=ΔmρF, color=colours[6], markersize=5, marker=:pentagon;label=label6)

            scatter!(plt3, [r], xerr = Δr, [w*mη],  yerr=Δproduct(mη,w,Δmη,Δw), color=colours[1] ,marker=:rect;label=label1)
            scatter!(plt3, [r], xerr = Δr, [w*ma],  yerr=Δproduct(ma,w,Δma,Δw), color=colours[2] ,marker=:rect;label=label2)
            scatter!(plt3, [r], xerr = Δr, [w*mπA], yerr=Δproduct(mπA,w,ΔmπA,Δw), color=colours[3], marker=:circ;label=label3)
            scatter!(plt3, [r], xerr = Δr, [w*mπF], yerr=Δproduct(mπF,w,ΔmπF,Δw), color=colours[4], marker=:circ;label=label4)
            scatter!(plt3, [r], xerr = Δr, [w*mρA], yerr=Δproduct(mρA,w,ΔmρA,Δw), color=colours[5] , marker=:pentagon;label=label5)
            scatter!(plt3, [r], xerr = Δr, [w*mρF], yerr=Δproduct(mρF,w,ΔmρF,Δw), color=colours[6] , marker=:pentagon;label=label6)

        end
        plot!(plt2,xlims=(0.87,0.92))
        plot!(plt3,xlims=(0.87,0.92))
    end

    savefig(plt0,joinpath(dir,"masses_vs_bare_mass.pdf"))
    savefig(plt1,joinpath(dir,"masses_vs_spatial_extent.pdf"))
    savefig(plt2,joinpath(dir,"masses_vs_mrho_over_mpi.pdf"))
    savefig(plt3,joinpath(dir,"masses_gradient_flow_scale.pdf"))

    return plt0, plt1, plt2
end