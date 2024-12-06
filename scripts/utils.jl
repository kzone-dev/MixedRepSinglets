eigenvalues_meff_mixed_rep(args...;kws...) = eigenvalues_eigenvectors_meff_mixed_rep(args...;kws...)[1:5]
function eigenvalues_eigenvectors_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    symmetry = +1 
    correlation_matrix = LatticeUtils._bin_correlator_matrix(correlation_matrix;binsize)
    if deriv 
        correlation_matrix = correlator_derivative(correlation_matrix;t_dim=4)
        symmetry = -1 
    end
    # use correlator binning
    eigenvalues_jackknife, eigenvectors_jackknife = eigenvalues_eigenvectors_jackknife_samples(correlation_matrix;t0)
    eigvals, Δeigvals = LatticeUtils.apply_jackknife(eigenvalues_jackknife;dims=2)
    eigvecs, Δeigvecs = LatticeUtils.apply_jackknife(eigenvectors_jackknife;dims=3)
    meff, Δmeff =  implicit_meff_jackknife(eigenvalues_jackknife;sign=symmetry)
    return eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife, eigvecs, Δeigvecs, eigenvectors_jackknife
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