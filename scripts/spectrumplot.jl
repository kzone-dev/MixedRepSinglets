function plot_spectrum(tablepath)
    results_MR = joinpath(tablepath,"table_results_MR.csv")

    # get meson masses
    data = readdlm(results_MR,';',skipstart=1)

    plt0 = plot()
    plt1 = plot()
    plt2 = plot()
    xticks0 = Float64[]
    xticks1 = Int[]

    for row in eachrow(data[:,3:18])
        T, L, mf, mas, ma, Δma, mη, Δmη, mπF, ΔmπF, mπA, ΔmπA, mρF, ΔmρF, mρA, ΔmρA = row

        r  = mπF/mρF
        Δr = sqrt((ΔmπF/mρF)^2 + (mπF*ΔmρF/mρF^2)^2)

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

    display(plt0)
    display(plt1)
    display(plt2)

    return plt0, plt1, plt2
end