function plot_spectrum(tablepath,plotdir,gradient_flow_results="input/gradient_flow_results.csv")
    results_MR = joinpath(tablepath,"table_results_MR.csv")

    dir = joinpath(plotdir,"spectrum")
    ispath(dir) || mkpath(dir)


    # get meson masses
    data = readdlm(results_MR,';',skipstart=1)
    gf = readdlm(gradient_flow_results,',',skipstart=1)

    plt0 = plot(legend=:right,ylabel=L"meson masses $[a^{-1}]$", xlabel=L"$m_0^{\rm f}$", title="meson masses vs. bare fermion mass")
    plt1 = plot(legend=:right,ylabel=L"meson masses $[a^{-1}]$", xlabel=L"N_t", title=L"meson masses vs. spatial extent $N_t$")
    plt2 = plot(legend=:right,ylabel=L"meson masses $[a^{-1}]$", xlabel=L"m_\pi^{\rm f} / m_\rho^{\rm f}", title = L"meson masses vs. fundamental $m_\pi / m_\rho ~ (N_t=64)$" )
    plt3 = plot(legend=:right,ylabel=L"meson masses $[w]$", xlabel=L"m_\pi^{\rm f} / m_\rho^{\rm f}", title=L"meson mass in gradient flow scale $w ~ (N_t=64)$")
    
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
        
        scatter!(plt0, [mf+offset], [ma], yerr=Δma,label="", marker=:rect, color=:red)
        scatter!(plt0, [mf+offset], [mη], yerr=Δmη,label="", marker=:rect, color=:red)
        scatter!(plt0, [mf+offset], [mπF], yerr=ΔmπF,label="", marker=:circ, color=:blue)
        scatter!(plt0, [mf+offset], [mπA], yerr=ΔmπA,label="", marker=:circ, color=:blue)
        scatter!(plt0, [mf+offset], [mρF], yerr=ΔmρF,label="", marker=:pentagon, color=:black)
        scatter!(plt0, [mf+offset], [mρA], yerr=ΔmρA,label="", marker=:pentagon, color=:black)
        plot!(plt0;xticks=xticks0)

        if isapprox(mf,-0.71)
            scatter!(plt1, [T], [ma], yerr=Δma,label="", marker=:rect, color=:red)
            scatter!(plt1, [T], [mη], yerr=Δmη,label="", marker=:rect, color=:red)
            scatter!(plt1, [T], [mπF], yerr=ΔmπF,label="", marker=:circ, color=:blue)
            scatter!(plt1, [T], [mπA], yerr=ΔmπA,label="", marker=:circ, color=:blue)
            scatter!(plt1, [T], [mρF], yerr=ΔmρF,label="", marker=:pentagon, color=:black)
            scatter!(plt1, [T], [mρA], yerr=ΔmρA,label="", marker=:pentagon, color=:black)
            plot!(plt1;xticks=xticks1)
        end

        if T == 64
            scatter!(plt2, [r], xerr = Δr, [ma], yerr=Δma,label="", marker=:rect, color=:red)
            scatter!(plt2, [r], xerr = Δr, [mη], yerr=Δmη,label="", marker=:rect, color=:red)
            scatter!(plt2, [r], xerr = Δr, [mπF], yerr=ΔmπF,label="", marker=:circ, color=:blue)
            scatter!(plt2, [r], xerr = Δr, [mπA], yerr=ΔmπA,label="", marker=:circ, color=:blue)
            scatter!(plt2, [r], xerr = Δr, [mρF], yerr=ΔmρF,label="", marker=:pentagon, color=:black)
            scatter!(plt2, [r], xerr = Δr, [mρA], yerr=ΔmρA,label="", marker=:pentagon, color=:black)

            scatter!(plt3, [r], xerr = Δr, [w*ma],  yerr=Δproduct(ma,w,Δma,Δw),label="", marker=:rect, color=:red)
            scatter!(plt3, [r], xerr = Δr, [w*mη],  yerr=Δproduct(mη,w,Δmη,Δw),label="", marker=:rect, color=:red)
            scatter!(plt3, [r], xerr = Δr, [w*mπF], yerr=Δproduct(mπF,w,ΔmπF,Δw),label="", marker=:circ, color=:blue)
            scatter!(plt3, [r], xerr = Δr, [w*mπA], yerr=Δproduct(mπA,w,ΔmπA,Δw),label="", marker=:circ, color=:blue)
            scatter!(plt3, [r], xerr = Δr, [w*mρF], yerr=Δproduct(mρF,w,ΔmρF,Δw),label="", marker=:pentagon, color=:black)
            scatter!(plt3, [r], xerr = Δr, [w*mρA], yerr=Δproduct(mρA,w,ΔmρA,Δw),label="", marker=:pentagon, color=:black)

        end
        plot!(plt2,xlims=(0.88,0.92))
        plot!(plt3,xlims=(0.88,0.92))
    end

    for plt in (plt0,plt1,plt2,plt3)
        scatter!(plt, [],[],label=L"$J^P = 0^-$ (singlet)",marker=:rect, color=:red)
        scatter!(plt, [],[],label=L"$J^P = 0^-$ (nonsinglet)",marker=:circ, color=:blue)
        scatter!(plt, [],[],label=L"$J^P = 1^+$ (nonsinglet)",marker=:pentagon, color=:black)
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