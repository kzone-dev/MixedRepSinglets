# Here, the uncertainties would naively lead to samples where you would try to take the FP64 square root of a negative number
# which leads to a domain error. By adding a vanishing imaginary part, julia will compute the copmplex result with a negligible 
# complex part, which I discard in the end-
function effective_mixing_angle(evecs_jackknife)
    # Note: The thirs entry denotes the i=1,2 entries of the eigenvectors, whereas the fourth index labels the eignvectors
    arg = @. - (evecs_jackknife[1,1,:,:] * evecs_jackknife[2,2,:,:]) / (evecs_jackknife[2,1,:,:] * evecs_jackknife[1,2,:,:])
    angle = atand.(sqrt.((arg .+ 0im) ))
    ϕ, Δϕ = MixedRepSinglets.apply_jackknife(angle,dims=1)
    N, T  = size(angle)
    ϕcov  = sqrt(N-1)*cov(real.(angle); dims=1,corrected=false)
    return real.(ϕ), real.(Δϕ), Hermitian(ϕcov)
end
function effective_mixing_angle_samples(evecs_jackknife)
    arg  = @. - (evecs_jackknife[1,1,:,:] * evecs_jackknife[2,2,:,:]) / (evecs_jackknife[2,1,:,:] * evecs_jackknife[1,2,:,:])
    ϕ_jk = atand.(sqrt.((arg) .+ 0im ))
    return real.(ϕ_jk) 
end
function _loss_of_signal(x,Δx;ind0=0,thresh)
    rel_error = @. abs(Δx[ind0+1:end]/x[ind0+1:end])
    ind = findfirst(x-> x > thresh, rel_error) + ind0
    return ind
end
function _plot_effective_mixing_angle(correlation_matrix,title;t0,binsize,deriv)
    T = last(size(correlation_matrix))
    evals, Δevals, meff, Δmeff, evals_jk, evecs, Δevecs, evecs_jk = eigenvalues_eigenvectors_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    ϕ, Δϕ, ϕcov =  effective_mixing_angle(evecs_jk)

    # indicate range of ground state signal
    t_meff = _loss_of_signal(meff[2,:],Δmeff[2,:];thresh=0.5)
    t_angle = _loss_of_signal(ϕ, Δϕ; ind0=t0, thresh=0.5)

    label ="effective mixing angle"
    label =""
    ylabel=L"$\phi/ ^\circ$ "
    xlabel=L"t > t_0 = %$t0"

    plt = plot(;title,ylabel,xlabel,ylims=(-5,20),xlims=(t0,T÷2))
    scatter!(plt,ϕ,yerr=Δϕ;label)
    vspan!(plt,[t_meff,T÷2],color=:grey,alpha=0.5,label="loss of signal in effective mass")
end
function _fit_effective_mixing_angle(correlation_matrix;t0,binsize,deriv)
    evals, Δevals, meff, Δmeff, evals_jk, evecs, Δevecs, evecs_jk = eigenvalues_eigenvectors_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    ϕ, Δϕ, ϕcov =  effective_mixing_angle(evecs_jk)
    # fit until we loose the signal in the effective mass of the ground state
    tmax = _loss_of_signal(meff[2,:],Δmeff[2,:];thresh=0.5) - 1
    tmin = t0 +1
    # extract fitted mxing angle 
    #ϕfit, Δϕfit = _fit_effective_mixing_angle(ϕ, Δϕ, tmin, tmax)
    ϕfit, Δϕfit = _fit_effective_mixing_angle(ϕ, ϕcov, tmin, tmax)
    @show ϕfit, Δϕfit
    return ϕfit, Δϕfit, tmin, tmax
end
function _fit_effective_mixing_angle_jackknife_error(correlation_matrix;t0,binsize,deriv,tmin=nothing,tmax=nothing)
    evals, Δevals, meff, Δmeff, evals_jk, evecs, Δevecs, evecs_jk = eigenvalues_eigenvectors_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    ϕ_jk  = effective_mixing_angle_samples(evecs_jk)
    ϕ, Δϕ, ϕcov =  effective_mixing_angle(evecs_jk)
    # I'm using π/2 periodicity of the mixing angle here
    # fit until we loose the signal in the effective mass of the ground state
    if isnothing(tmax) 
        tmax = _loss_of_signal(meff[2,:],Δmeff[2,:];thresh=0.5) - 1
    end
    if isnothing(tmin)
        tmin = t0 +1
    end

    # fit individual jackknife samples
    # always use the jackknife estimate for the uncertainty as the weight in fitting
    N  = first(size(ϕ_jk))
    ϕs = zeros(N) 
    for i in 1:N
        #ϕs[i] = first(_fit_effective_mixing_angle(ϕ_jk[i,:], Δϕ, tmin, tmax))
        ϕs[i] = first(_fit_effective_mixing_angle(ϕ_jk[i,:], ϕcov, tmin, tmax))
    end
    ϕfit, Δϕfit = MixedRepSinglets.apply_jackknife(ϕs)
    return ϕfit, Δϕfit, tmin, tmax
end
function _fit_effective_mixing_angle(ϕ, Δϕ::AbstractVector, tmin, tmax)
    # fit to a constant
    @. model(x,p) = p[1] + 0*x
    p0  = [0.0] # initial guess: no mixing 
    t   = tmin:tmax
    fit = curve_fit(model,t,ϕ[t],(Δϕ[t]).^(-2),p0) 
    # extract fitted mxing angle 
    ϕfit, Δϕfit = fit.param[1], stderror(fit)[1]
    return ϕfit, Δϕfit
end
function _fit_effective_mixing_angle(ϕ, ϕcov::AbstractMatrix, tmin, tmax)
    # fit to a constant
    @. model(x,p) = p[1] + 0*x
    p0  = [0.0] # initial guess: no mixing 
    t   = tmin:tmax
    fit = curve_fit(model,t,ϕ[t],Hermitian(inv(ϕcov[t,t])),p0) 
    # extract fitted mxing angle 
    ϕfit, Δϕfit = fit.param[1], stderror(fit)[1]
    return ϕfit, Δϕfit
end
function _plot_title(h5corrs,ensemble)
    β   = h5read(h5corrs,joinpath(ensemble,"beta"))
    T,L = h5read(h5corrs,joinpath(ensemble,"lattice"))[1:2]
    mf  = h5read(h5corrs,joinpath(ensemble,"quarkmasses_fundamental"))[1]
    mas = h5read(h5corrs,joinpath(ensemble,"quarkmasses_antisymmetric"))[1]
    title = L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"
    return title
end
function plot_and_write_mixing_angles(parameters_gevp,hdf5path,tablepath,tex_tablepath,plotpath)
    h5corrs   = joinpath(hdf5path,"singlets_smeared_correlators.hdf5")

    file_mixing = joinpath(tablepath,"table_mixing_angle.csv")
    file_mixing_MR = joinpath(tablepath,"table_mixing_angle_MR.csv")
    io_mixing = open(file_mixing,"w")
    io_mixing_MR = open(file_mixing_MR,"w")

    write(io_mixing_MR,L"Label;$\beta$;$N_t$;$N_l$;$\phi$;$\Delta \phi$","\n")
    write(io_mixing,L"Label;$~~~~\beta~~~~$;$~~~~N_t~~~~$;$~~~~N_s~~~~$;$~~~~\phi/{}^{\circ}~~~~$","\n")
    plotpath = joinpath(plotpath,"mixing_angle")

    parameters = readdlm(parameters_gevp,';';skipstart=1)
    for row in eachrow(parameters)

        ensemble, channel, t0, binsize, = row[1:4]
        nops, deriv  = [1,10], false
        t0  = 5

        channel == "g5_singlet" || continue    
        matrixname ="correlation_matrix_g5_singlet"
        correlation_matrix = h5read(h5corrs,joinpath(ensemble,matrixname))
        correlation_matrix = correlation_matrix[nops,nops,:,:]

        β   = h5read(h5corrs,joinpath(ensemble,"beta"))
        T,L = h5read(h5corrs,joinpath(ensemble,"lattice"))[1:2]
        mf  = h5read(h5corrs,joinpath(ensemble,"quarkmasses_fundamental"))[1]
        mas = h5read(h5corrs,joinpath(ensemble,"quarkmasses_antisymmetric"))[1]

        # special case M1, where the correlated fitr is not stable
        tmax = ensemble == "M1" ? 11 : nothing
        # get fitted mixing angle 
        ϕ, Δϕ, t1, t2 = _fit_effective_mixing_angle_jackknife_error(correlation_matrix;t0,binsize,deriv,tmax)
        title = "ensemble $ensemble" #_plot_title(h5corrs,ensemble)    
        plt   = _plot_effective_mixing_angle(correlation_matrix,title;t0,binsize,deriv)
        # add best fit to effective mixing angle
        add_fit_range!(plt,t1,t2,ϕ,Δϕ;label=L"fit: mixing angle $\phi/ ^\circ$")

        write(io_mixing_MR,"$ensemble;$β;$T;$L;$ϕ;$Δϕ","\n")
        write(io_mixing,"$ensemble;$β;$T;$L;$(errorstring(ϕ,Δϕ;nsig=2))","\n")

        display(plt)
        ispath(plotpath) || mkpath(plotpath)
        savefig(plt, joinpath(plotpath,"mixing_angle_$ensemble.pdf"))
    end
    close(io_mixing)
    close(io_mixing_MR)

    # convert table into a tex compatible formatting
    mixing = readdlm(joinpath(tablepath,"table_mixing_angle.csv"),';')
    write_tex_table(joinpath(tex_tablepath,"table_mixing.tex"),mixing)
end
