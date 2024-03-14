function plot_spectrum(tablepath)
    results_MR = joinpath(tablepath,"table_results_MR.csv")

    # get meson masses
    data = readdlm(results_MR,';',skipstart=1)

    plt0 = plot(legend=:right,ylabel=L"meson masses $[a^{-1}]$", xlabel=L"$m_0^{\rm f}$", title="meson masses vs. bare fermion mass")
    plt1 = plot(legend=:right,ylabel=L"meson masses $[a^{-1}]$", xlabel=L"N_t", title=L"meson masses vs. spatial extent $N_t$")
    plt2 = plot(legend=:right,ylabel=L"meson masses $[a^{-1}]$", xlabel=L"m_\pi / m_\rho", title = L"meson masses vs. $m_\pi / m_\rho ~ (N_t=64)$" )
    #plt3 = plot(legend=:right,ylabel=L"meson masses $[m_\pi^{(f)}]$", xlabel=L"m_\pi / m_\rho", title=L"meson mass ratios $(N_t=64)$")
    
    xticks0 = Float64[]
    xticks1 = Int[]

    for row in eachrow(data[:,3:18])
        T, L, mf, mas, ma, Δma, mη, Δmη, mπF, ΔmπF, mπA, ΔmπA, mρF, ΔmρF, mρA, ΔmρA = row

        Δratio(x,y,Δx,Δy) = sqrt((Δx/y)^2 + (x*Δy/y^2)^2) 
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
        end
        plot!(plt2,xlims=(0.88,0.92))
    end

    for plt in (plt1,plt2,plt0)
        scatter!(plt, [],[],label=L"$J^P = 0^-$ (singlet)",marker=:rect, color=:red)
        scatter!(plt, [],[],label=L"$J^P = 0^-$ (nonsinglet)",marker=:circ, color=:blue)
        scatter!(plt, [],[],label=L"$J^P = 1^+$ (nonsinglet)",marker=:pentagon, color=:black)
    end

    display(plt0)
    display(plt1)
    display(plt2)
    return plt0, plt1, plt2
end