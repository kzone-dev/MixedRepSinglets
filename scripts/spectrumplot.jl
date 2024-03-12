using DelimitedFiles
using Plots
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=11, tickfontsize=10, labelfontsize=14, markersize=5)

# get meson masses
data = readdlm("output/tables/table_results_MR.csv",';',skipstart=1)

plt0 = plot()
plt1 = plot()
plt2 = plot()
xticks = Float64[]

for row in eachrow(data[:,3:18])
    T, L, mf, mas, ma, Δma, mη, Δmη, mπF, ΔmπF, mπA, ΔmπA, mρF, ΔmρF, mρA, ΔmρA = row

    r  = mπF/mρF
    Δr = ΔmπF/mρF + mπF*ΔmρF/mρF^2

    offset = sign(T-64)*0.0005
    push!(xticks,mf)
    unique!(xticks)
    
    scatter!(plt0, [mf+offset], [ma], yerr=Δma,label="", marker=:rect, color=:red)
    scatter!(plt0, [mf+offset], [mη], yerr=Δmη,label="", marker=:rect, color=:red)
    scatter!(plt0, [mf+offset], [mπF], yerr=ΔmπF,label="", marker=:circ, color=:blue)
    scatter!(plt0, [mf+offset], [mπA], yerr=ΔmπA,label="", marker=:circ, color=:blue)
    scatter!(plt0, [mf+offset], [mρF], yerr=ΔmρF,label="", marker=:pentagon, color=:black)
    scatter!(plt0, [mf+offset], [mρA], yerr=ΔmρA,label="", marker=:pentagon, color=:black)
    plot!(plt0;xticks)

    if isapprox(mf,-0.71)
        scatter!(plt1, [T], [ma], yerr=Δma,label="", marker=:rect, color=:red)
        scatter!(plt1, [T], [mη], yerr=Δmη,label="", marker=:rect, color=:red)
        scatter!(plt1, [T], [mπF], yerr=ΔmπF,label="", marker=:circ, color=:blue)
        scatter!(plt1, [T], [mπA], yerr=ΔmπA,label="", marker=:circ, color=:blue)
        scatter!(plt1, [T], [mρF], yerr=ΔmρF,label="", marker=:pentagon, color=:black)
        scatter!(plt1, [T], [mρA], yerr=ΔmρA,label="", marker=:pentagon, color=:black)
    end

    scatter!(plt2, [r], xerr = Δr, [ma], yerr=Δma,label="", marker=:rect, color=:red)
    scatter!(plt2, [r], xerr = Δr, [mη], yerr=Δmη,label="", marker=:rect, color=:red)
    scatter!(plt2, [r], xerr = Δr, [mπF], yerr=ΔmπF,label="", marker=:circ, color=:blue)
    scatter!(plt2, [r], xerr = Δr, [mπA], yerr=ΔmπA,label="", marker=:circ, color=:blue)
    scatter!(plt2, [r], xerr = Δr, [mρF], yerr=ΔmρF,label="", marker=:pentagon, color=:black)
    scatter!(plt2, [r], xerr = Δr, [mρA], yerr=ΔmρA,label="", marker=:pentagon, color=:black)

end

plt0
plt1
plt2