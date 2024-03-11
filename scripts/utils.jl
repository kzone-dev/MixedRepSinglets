function _copy_lattice_parameters(outfile,infile,ensemble)
    file = h5open(infile)[ensemble]
    entries = filter(!contains("correlation_matrix"),keys(file))
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
function write_eigenvalues_and_effective_masses(correlation_matrix,outputfile,inputfile,ensemble,name; t0 = 1, binsize = 2, deriv = true, setup = true)
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    setup && _copy_lattice_parameters(outputfile,inputfile,ensemble)

    new_name = replace(name,"correlation_matrix_" => "")
    h5write(outputfile,joinpath(ensemble,"meff_$new_name"),meff)
    h5write(outputfile,joinpath(ensemble,"eigvals_$new_name"),eigvals)
    h5write(outputfile,joinpath(ensemble,"Delta_meff_$new_name"),Δmeff)
    h5write(outputfile,joinpath(ensemble,"Delta_eigvals_$new_name"),Δeigvals)
    # generic quantitites used in the GEVP inversion and data preparation
    h5write(outputfile,joinpath(ensemble,"t0_$new_name"),t0)
    h5write(outputfile,joinpath(ensemble,"deriv_$new_name"),deriv)
    h5write(outputfile,joinpath(ensemble,"binsize_$new_name"),binsize)
end