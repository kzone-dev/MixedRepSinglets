function _copy_lattice_parameters(outfile,infile,ensemble;group="")
    file = h5open(infile)[ensemble]
    entries = filter(!contains("correlation_matrix"),keys(file))
    for entry in entries
        label = joinpath(ensemble,group,entry)
        h5write(outfile,label,read(file,entry))
    end
end
function eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    symmetry = +1 
    correlation_matrix = correlator_folding(correlation_matrix;t_dim=4,sign=symmetry)
    correlation_matrix = _bin_correlator_matrix(correlation_matrix;binsize)
    #correlation_matrix = correlation_matrix[:,:,1:binsize:end,:]
    if deriv 
        correlation_matrix = correlator_derivative(correlation_matrix;t_dim=4)
        symmetry = -1 
    end
    # use correlator binning
    eigvals, Δeigvals = eigenvalues(correlation_matrix;t0)
    eigenvalues_jackknife = eigenvalues_jackknife_samples(correlation_matrix;t0)
    meff, Δmeff =  meff_from_jackknife(eigenvalues_jackknife;sign=symmetry,swap=nothing)
    return eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife
end
function write_eigenvalues_and_effective_masses(correlation_matrix,outputfile,inputfile,ensemble,channel; t0, binsize, deriv, resamples = false)
    eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
   
    _copy_lattice_parameters(outputfile,inputfile,ensemble;group=channel)

    h5write(outputfile,joinpath(ensemble,channel,"meff"),meff)
    h5write(outputfile,joinpath(ensemble,channel,"eigvals"),eigvals)
    h5write(outputfile,joinpath(ensemble,channel,"Delta_meff"),Δmeff)
    h5write(outputfile,joinpath(ensemble,channel,"Delta_eigvals"),Δeigvals)

    # generic quantitites used in the GEVP inversion and data preparation
    h5write(outputfile,joinpath(ensemble,channel,"t0"),t0)
    h5write(outputfile,joinpath(ensemble,channel,"deriv"),deriv)
    h5write(outputfile,joinpath(ensemble,channel,"binsize"),binsize)

    if resamples
        h5write(outputfile,joinpath(ensemble,channel,"eigvals_resamples"),eigenvalues_jackknife)
        h5write(outputfile,joinpath(ensemble,channel,"eigvals_resample_type"),"jackknife")
    end

end
function _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,β,T,L,mf,mas;nstates=1,tmax=nothing,tag="")
    Nops = first(size(meff))

    title = L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"   
    plt1 = plot(title=title, xlabel="t", ylabel="effective mass")
    plt2 = plot(title=title, xlabel="t", ylabel="eigenvalues")

    function taggedlabel(state)
        label = isequal(state,0) ? "ground state" : "excited state #$(state)" 
        label = isempty(tag) ? label : "$tag: $label"
        return label
    end
    
    for state in 0:nstates-1
        if isnothing(tmax) 
            rel_error = abs.(Δmeff[Nops-state,:]./meff[Nops-state,:])
            tmax = findfirst(x-> x > (1/2), rel_error)
            tmax = isnothing(tmax) ? T÷2 : tmax
        end
        range = 2:tmax
        scatter!(plt1,range, meff[Nops-state,range], yerr= Δmeff[Nops-state,range],label=taggedlabel(state))
        plot_correlator!(plt2,range,eigvals[Nops-state,range],Δeigvals[Nops-state,range];label=taggedlabel(state))
    end
    return plt1, plt2
end