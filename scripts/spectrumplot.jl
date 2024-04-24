function plot_spectrum(tablepath,plotdir,gradient_flow_results="input/gradient_flow_results.csv")
    results_MR = joinpath(tablepath,"table_results_MR.csv")

    dir = joinpath(plotdir,"spectrum")
    ispath(dir) || mkpath(dir)

    colours = palette(:Paired_6)

    # get meson masses
    data = readdlm(results_MR,';',skipstart=1)
    gf = readdlm(gradient_flow_results,',',skipstart=1)

    plt0 = plot(legend_column=2,legendfontsize=13,legend=:left,ylabel=L"meson masses $[a^{-1}]$", xlabel=L"$m_0^{\rm f}$", title="meson masses vs. bare fermion mass")
    plt1 = plot(legend_column=2,legendfontsize=13,legend=:left,ylabel=L"meson masses $[a^{-1}]$", xlabel=L"N_t", title=L"meson masses vs. spatial extent $N_t$")
    plt2 = plot(legend_column=2,legendfontsize=13,legend=:left,ylabel=L"meson masses $[a^{-1}]$", xlabel=L"m_\pi^{\rm f} / m_\rho^{\rm f}", title = L"meson masses vs. fundamental $m_\pi / m_\rho ~ (N_t=64)$" )
    plt3 = plot(legend_column=2,legendfontsize=13,legend=:left,ylabel=L"meson masses $[w]$", xlabel=L"m_\pi^{\rm f} / m_\rho^{\rm f}", title=L"meson mass in gradient flow scale $w ~ (N_t=64)$")
    
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
        
        scatter!(plt0, [mf+offset], [mη], yerr=Δmη,label="", markersize=5, marker=:rect, color=colours[1])
        scatter!(plt0, [mf+offset], [ma], yerr=Δma,label="", markersize=5, marker=:rect, color=colours[2])
        scatter!(plt0, [mf+offset], [mπA], yerr=ΔmπA,label="", markersize=5, marker=:circ, color=colours[3])
        scatter!(plt0, [mf+offset], [mπF], yerr=ΔmπF,label="", markersize=5, marker=:circ, color=colours[4])
        scatter!(plt0, [mf+offset], [mρA], yerr=ΔmρA,label="", markersize=5, marker=:pentagon, color=colours[5])
        scatter!(plt0, [mf+offset], [mρF], yerr=ΔmρF,label="", markersize=5, marker=:pentagon, color=colours[6])
        plot!(plt0;xticks=xticks0)

        if isapprox(mf,-0.71)
            scatter!(plt1, [T], [mη], yerr=Δmη, color=colours[1], label="",markersize=5, marker=:rect)
            scatter!(plt1, [T], [ma], yerr=Δma, color=colours[2], label="",markersize=5, marker=:rect)
            scatter!(plt1, [T], [mπA], yerr=ΔmπA, color=colours[3], label="",markersize=5, marker=:circ)
            scatter!(plt1, [T], [mπF], yerr=ΔmπF, color=colours[4], label="",markersize=5, marker=:circ)
            scatter!(plt1, [T], [mρA], yerr=ΔmρA, color=colours[5], label="",markersize=5, marker=:pentagon)
            scatter!(plt1, [T], [mρF], yerr=ΔmρF, color=colours[6], label="",markersize=5, marker=:pentagon)
            plot!(plt1;xticks=xticks1)
        end

        if T == 64
            scatter!(plt2, [r], xerr = Δr, [mη], yerr=Δmη, color=colours[1] ,label="", markersize=5, marker=:rect)
            scatter!(plt2, [r], xerr = Δr, [ma], yerr=Δma, color=colours[2] ,label="", markersize=5, marker=:rect)
            scatter!(plt2, [r], xerr = Δr, [mπA], yerr=ΔmπA, color=colours[3] ,label="", markersize=5, marker=:circ)
            scatter!(plt2, [r], xerr = Δr, [mπF], yerr=ΔmπF, color=colours[4] ,label="", markersize=5, marker=:circ)
            scatter!(plt2, [r], xerr = Δr, [mρA], yerr=ΔmρA, color=colours[5] ,label="", markersize=5, marker=:pentagon)
            scatter!(plt2, [r], xerr = Δr, [mρF], yerr=ΔmρF, color=colours[6] ,label="", markersize=5, marker=:pentagon)

            scatter!(plt3, [r], xerr = Δr, [w*mη],  yerr=Δproduct(mη,w,Δmη,Δw), color=colours[1] ,label="", marker=:rect)
            scatter!(plt3, [r], xerr = Δr, [w*ma],  yerr=Δproduct(ma,w,Δma,Δw), color=colours[2] ,label="", marker=:rect)
            scatter!(plt3, [r], xerr = Δr, [w*mπA], yerr=Δproduct(mπA,w,ΔmπA,Δw), color=colours[3] ,label="", marker=:circ)
            scatter!(plt3, [r], xerr = Δr, [w*mπF], yerr=Δproduct(mπF,w,ΔmπF,Δw), color=colours[4] ,label="", marker=:circ)
            scatter!(plt3, [r], xerr = Δr, [w*mρA], yerr=Δproduct(mρA,w,ΔmρA,Δw), color=colours[5] ,label="", marker=:pentagon)
            scatter!(plt3, [r], xerr = Δr, [w*mρF], yerr=Δproduct(mρF,w,ΔmρF,Δw), color=colours[6] ,label="", marker=:pentagon)

        end
        plot!(plt2,xlims=(0.87,0.92))
        plot!(plt3,xlims=(0.87,0.92))
    end

    for plt in (plt0,plt1,plt2,plt3)
        scatter!(plt, [],[],label=L"\eta'",marker=:rect, color=colours[1])
        scatter!(plt, [],[],label=L"a",marker=:rect, color=colours[2])
        scatter!(plt, [],[],label=L"ps",marker=:circ, color=colours[4])
        scatter!(plt, [],[],label=L"PS",marker=:circ, color=colours[3])
        scatter!(plt, [],[],label=L"v",marker=:pentagon, color=colours[6])
        scatter!(plt, [],[],label=L"V",marker=:pentagon, color=colours[5])
    end

    display(plt0)
    display(plt1)
    display(plt2)
    display(plt3)
    savefig(plt0,joinpath(dir,"masses_vs_bare_mass.pdf"))
    savefig(plt1,joinpath(dir,"masses_vs_spatial_extent.pdf"))
    savefig(plt2,joinpath(dir,"masses_vs_mrho_over_mpi.pdf"))
    savefig(plt3,joinpath(dir,"masses_gradient_flow_scale.pdf"))

    return plt0, plt1, plt2
end