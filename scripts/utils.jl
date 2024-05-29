function run_corrfitter(parameters,hdf5file,outdir)
    args = `$(abspath(parameters)) $(abspath(hdf5file)) $(abspath(outdir))`
    try
        run(`python3 $(abspath("scripts/fitting_eigenvalues.py")) $args`)
    catch
        run(`python  $(abspath("scripts/fitting_eigenvalues.py")) $args`)
    end
end
function main_write_hdf5_logs(path,h5file,parameterfile;regexp=false)
    input = readdlm(parameterfile,';',skipstart=1)
    for prm in eachrow(input)
        
        dir, file, typeCONN, rep, name = prm
        fileCONN = joinpath(path,dir,"out",file)
        @show fileCONN   

        if regexp 
            typeCONN = Regex(typeCONN)
            writehdf5_spectrum_with_regexp(fileCONN,h5file,typeCONN,h5group="$name/$rep/CONN",sort=false)
        else
            writehdf5_spectrum(fileCONN,h5file,typeCONN,h5group="$name/$rep/CONN",sort=false)
        end
        
    end    
end
function eigenvalues_meff(correlation_matrix;t0=1,binsize=1,deriv=false)
    symmetry = +1 
    correlation_matrix = correlator_folding(correlation_matrix;t_dim=4,sign=symmetry)
    correlation_matrix = _bin_correlator_matrix(correlation_matrix;binsize)

    if deriv 
        correlation_matrix = correlator_derivative(correlation_matrix;t_dim=4)
        symmetry = -1 
    end

    eigvals, Δeigvals = eigenvalues(correlation_matrix;t0)
    eigenvalues_jackknife = eigenvalues_jackknife_samples(correlation_matrix;t0)
    meff, Δmeff =  meff_from_jackknife(eigenvalues_jackknife;sign=symmetry,swap=nothing)
    return eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife
end
