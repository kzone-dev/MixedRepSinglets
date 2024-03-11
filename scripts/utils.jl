function _copy_lattice_parameters(outfile,infile,ensemble)
    file = h5open(infile)[ensemble]
    entries = filter(!contains("correlation_matrix_g5_singlet"),keys(file))
    for entry in entries
        h5write(outfile,joinpath(ensemble,entry),read(file,entry))
    end
end
function eigenvalues_meff_mixed_rep(correlation_matrix;t0 = 1, binsize = 1, deriv = true)
    symmetry = +1 
    correlation_matrix = correlator_folding(correlation_matrix;t_dim=4,sign=symmetry)
    correlation_matrix = _bin_correlator_matrix(correlation_matrix;binsize)
    if deriv 
        correlation_matrix = correlator_derivative(correlation_matrix;t_dim=4)
        symmetry = -1 
    end
    # use correlator binning
    eigvals, Δeigvals = eigenvalues(correlation_matrix;t0)
    eigenvalues_jackknife = eigenvalues_jackknife_samples(correlation_matrix;t0 ,imag_thresh = 1E-12)
    meff, Δmeff =  meff_from_jackknife(eigenvalues_jackknife;sign=symmetry,swap=nothing)
    return eigvals, Δeigvals, meff, Δmeff
end
function write_eigenvalues_and_effective_masses(outputfile,inputfile,ensemble,name; t0 = 1, binsize = 2, deriv = true)
    correlation_matrix = h5read(h5corrs,joinpath(ensemble,name))
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    _copy_lattice_parameters(outputfile,inputfile,ensemble)
    h5write(outputfile,joinpath(ensemble,"eigvals_g5_singlet"),eigvals)
    h5write(outputfile,joinpath(ensemble,"Delta_eigvals_g5_singlet"),Δeigvals)
    h5write(outputfile,joinpath(ensemble,"meff_g5_singlet"),meff)
    h5write(outputfile,joinpath(ensemble,"Delta_meff_g5_singlet"),Δmeff)
    h5write(outputfile,joinpath(ensemble,"t0"),t0)
    h5write(outputfile,joinpath(ensemble,"deriv"),deriv)
    h5write(outputfile,joinpath(ensemble,"binsize"),binsize)
end