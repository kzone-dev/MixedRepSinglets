using Pkg; Pkg.activate(".")
using MixedRepSinglets
using LaTeXStrings
using DelimitedFiles
using HDF5
using Statistics
using MixedRepSinglets
using Plots
gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=12, tickfontsize=12, labelfontsize=14, titlefontsize=14,  markersize=5)
include("utils.jl")
include("eigenvalues.jl")


function effective_mixing_angle(evecs_jackknife)
    arg = @. evecs_jackknife[2,1,:,:] * evecs_jackknife[1,2,:,:] / evecs_jackknife[1,1,:,:] / evecs_jackknife[2,2,:,:]
    angle = atand.(sqrt.(abs.(arg) ))
    ϕ, Δϕ = MixedRepSinglets.apply_jackknife(angle,dims=1)
    return ϕ, Δϕ 
end
function _loss_of_signal(x,Δx;ind0=0,thresh)
    rel_error = @. abs(Δx[ind0+1:end]/x[ind0+1:end])
    @show size(rel_error)
    ind = findfirst(x-> x > thresh, rel_error) + ind0
    return ind
end
function _plot_effective_mixing_angle(correlation_matrix,title;t0,binsize,deriv)
    T = last(size(correlation_matrix))
    evals, Δevals, meff, Δmeff, evals_jk, evecs, Δevecs, evecs_jk = eigenvalues_eigenvectors_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    ϕ, Δϕ =  effective_mixing_angle(evecs_jk)
    # I'm using π/2 periodicity of the mixing angle here
    ϕ = ϕ .- 90

    # indicate range of ground state signal
    t_meff = _loss_of_signal(meff[2,:],Δmeff[2,:];thresh=0.5)
    t_angle = _loss_of_signal(ϕ, Δϕ; ind0=t0, thresh=0.5)
    @show t_angle

    label ="effective mixing angle"
    label =""
    ylabel=L"effective mixing angle $\phi [°]$ "
    xlabel=L"t > t_0 = %$t0"

    plt = plot(;title,ylabel,xlabel,ylims=(-15,15),xlims=(t0,T÷2))
    scatter!(plt,ϕ,yerr=Δϕ;label)
    vspan!(plt,[t_meff,T÷2],color=:grey,alpha=0.5,label="loss of signal in effective mass")
end
function _plot_title(h5corrs,ensemble)
    β   = h5read(h5corrs,joinpath(ensemble,"beta"))
    T,L = h5read(h5corrs,joinpath(ensemble,"lattice"))[1:2]
    mf  = h5read(h5corrs,joinpath(ensemble,"quarkmasses_fundamental"))[1]
    mas = h5read(h5corrs,joinpath(ensemble,"quarkmasses_antisymmetric"))[1]
    title = L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"
    return title
end

parameters_gevp = joinpath("input","parameters_gevp.csv")
hdf5path    = "/home/fabian/Downloads/hdf5out_modified"
h5corrs     = joinpath(hdf5path,"singlets_smeared_correlators.hdf5")
h5eigenvals = joinpath(hdf5path,"singlets_smeared_eigenvalues.hdf5")

parameters = readdlm(parameters_gevp,';';skipstart=1)
for row in eachrow(parameters)

    ensemble, channel, t0, binsize, = row[1:4]
    nops, deriv  = [1,10], false
    t0  = 5

    channel == "g5_singlet" || continue    
    matrixname ="correlation_matrix_g5_singlet"
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,matrixname))
    correlation_matrix = correlation_matrix[nops,nops,:,:]

    title = _plot_title(h5corrs,ensemble)    
    plt   = _plot_effective_mixing_angle(correlation_matrix,title;t0,binsize,deriv)
    display(plt)
    plotpath = joinpath("plot","mixing_angle")
    ispath(plotpath) || mkpath(plotpath)
    savefig(plt, joinpath(plotpath,"mixing_angle_$ensemble.pdf"))
end
