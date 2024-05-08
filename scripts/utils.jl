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
function eigenvalues_eigenvectors_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    symmetry = +1 
    correlation_matrix = correlator_folding(correlation_matrix;t_dim=4,sign=symmetry)
    #correlation_matrix = _bin_correlator_matrix(correlation_matrix;binsize)
    correlation_matrix = correlation_matrix[:,:,1:binsize:end,:]
    if deriv 
        correlation_matrix = correlator_derivative(correlation_matrix;t_dim=4)
        symmetry = -1 
    end
    # use correlator binning
    eigvals, Δeigvals, eigvecs, Δeigvecs = eigenvalues_eigenvectors(correlation_matrix;t0)
    eigenvalues_jackknife, eigenvectors_jackknife = eigenvalues_eigenvectors_jackknife_samples(correlation_matrix;t0)
    meff, Δmeff =  meff_from_jackknife(eigenvalues_jackknife;sign=symmetry,swap=nothing)
    return eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife, eigvecs, Δeigvecs, eigenvectors_jackknife
end
function write_eigenvalues_and_effective_masses(correlation_matrix,outputfile,inputfile,ensemble,channel; t0, binsize, deriv, resamples = false)
    eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    eigvals_cov = MixedRepSinglets.cov_jackknife_eigenvalues(eigenvalues_jackknife)
    
    _copy_lattice_parameters(outputfile,inputfile,ensemble;group=channel)

    h5write(outputfile,joinpath(ensemble,channel,"meff"),meff)
    h5write(outputfile,joinpath(ensemble,channel,"eigvals"),eigvals)
    h5write(outputfile,joinpath(ensemble,channel,"Delta_meff"),Δmeff)
    h5write(outputfile,joinpath(ensemble,channel,"Delta_eigvals"),Δeigvals)
    h5write(outputfile,joinpath(ensemble,channel,"eigvals_cov"),eigvals_cov)

    # generic quantitites used in the GEVP inversion and data preparation
    h5write(outputfile,joinpath(ensemble,channel,"t0"),t0)
    h5write(outputfile,joinpath(ensemble,channel,"deriv"),deriv)
    h5write(outputfile,joinpath(ensemble,channel,"binsize"),binsize)

    if resamples
        h5write(outputfile,joinpath(ensemble,channel,"eigvals_resamples"),eigenvalues_jackknife)
        h5write(outputfile,joinpath(ensemble,channel,"eigvals_resample_type"),"jackknife")
    end

end
function channel_label(channel,state)    
    tag ="" # empty fallback
    
    if isequal(channel,"g5_singlet")
        #state == 0 && return L"a" 
        #state == 1 && return L"\eta'"
        state == 0 && return L"\eta^{\prime}_l" 
        state == 1 && return L"\eta^{\prime}_h"
    end

    isequal(channel,"g5_nonsinglet_FUN") && (tag = L"{\rm PS}")
    isequal(channel,"g5_nonsinglet_AS")  && (tag = L"{\rm ps}")
    isequal(channel,"g1_nonsinglet_FUN") && (tag = L"{\rm V}")
    isequal(channel,"g1_nonsinglet_AS")  && (tag = L"{\rm v}")

    label = isequal(state,0) ? "ground state" : "excited state #$(state)" 
    return L"%$tag (%$label)"
end
function _plot_meff_eigvals(meff,Δmeff,eigvals,Δeigvals,channel;nstates=1,tmax=nothing,title="",kws...)
    Nops, T = size(meff)

    plt1 = plot(title=title, xlabel="t", ylabel=L"am^{\rm eff}")
    plt2 = plot(title=title, xlabel="t", ylabel="eigenvalues")

    markershapes = (:circle, :diamond, :rect, :pentagon,  :octagon)
    
    for state in 0:nstates-1
        if isnothing(tmax) 
            thresh = !contains(channel,"nonsinglet") ? 0.5 : 0.25
            rel_error = abs.(Δmeff[Nops-state,:]./meff[Nops-state,:])
            tmax = findfirst(x-> x > thresh, rel_error)
            tmax = isnothing(tmax) ? T÷2 : tmax
        end
        ms = markershapes[state+1]
        range = 2:tmax
        scatter!(plt1,range, meff[Nops-state,range], yerr= Δmeff[Nops-state,range],label=channel_label(channel,state),markershape=ms,kws...)
        plot_correlator!(plt2,1:T,eigvals[Nops-state,1:T],Δeigvals[Nops-state,1:T];label=channel_label(channel,state),markershape=ms,kws...)
        tmax = nothing
    end
    return plt1, plt2
end
function errorstring(x,Δx;nsig=2)
    @assert Δx > 0
    sgn = x < 0 ? "-" : ""
    x = abs(x)  
    # round error part to desired number of signficant digits
    # convert to integer if no fractional part exists
    Δx_rounded = round(Δx,sigdigits=nsig) 
    # get number of decimal digits for x  
    floor_log_10 = floor(Int,log10(Δx))
    dec_digits   = (nsig - 1) - floor_log_10
    # round x, to desired number of decimal digits 
    # (standard julia function deals with negative dec_digits) 
    x_rounded = round(x,digits=dec_digits)
    # get decimal and integer part if there is a decimal part
    if dec_digits > 0
        digits_val = Int(round(x_rounded*10.0^(dec_digits)))
        digits_unc = Int(round(Δx_rounded*10.0^(dec_digits)))
        str_val = _insert_decimal(digits_val,dec_digits) 
        str_unc = _insert_decimal(digits_unc,dec_digits)
        str_unc = nsig > dec_digits ? str_unc : string(digits_unc)
        return sgn*"$str_val($str_unc)"
    else
        return sgn*"$(Int(x_rounded))($(Int(Δx_rounded)))"
    end
end
function _insert_decimal(val::Int,digits)
    str = lpad(string(val),digits,"0")
    pos = length(str) - digits
    int = rpad(str[1:pos],1,"0")
    dec = str[pos+1:end]
    inserted = int*"."*dec
    return inserted
end