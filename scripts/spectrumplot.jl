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