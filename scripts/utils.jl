function run_corrfitter(parameters,hdf5file,outdir;fpi=false)
    if fpi
        args = `$(abspath(parameters)) $(abspath(hdf5file)) $(abspath(outdir)) fpi`
    else    
        args = `$(abspath(parameters)) $(abspath(hdf5file)) $(abspath(outdir)) normal`
    end
    # execute python script
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
function _effective_mass_plots(results,fitparam)
    for (i,line) in enumerate(eachrow(results))
        ens, channel, rep, T, L, β, m, Δm, χ2dof = line
        tmin, tmax, tp, Nmax = fitparam[i,4:7]

        title = "$T × $L^3, β=$β, $channel, $rep"
        
        if channel == "g1"
            label1 = "$ens/$rep/CONN/DEFAULT_SEMWALL TRIPLET/g1"
            label2 = "$ens/$rep/CONN/DEFAULT_SEMWALL TRIPLET/g2"
            label3 = "$ens/$rep/CONN/DEFAULT_SEMWALL TRIPLET/g3"
            corr = (h5read(h5file,label1) .+ h5read(h5file,label2) .+ h5read(h5file,label3)) ./ 3 
        else
            label = "$ens/$rep/CONN/DEFAULT_SEMWALL TRIPLET/$channel"
            corr = h5read(h5file,label)
        end

        corr = correlator_folding(corr;t_dim=2,sign=+1)
        meff, Δmeff = implicit_meff_jackknife(corr')

        range = 1:div(T,2)
        ylims = (m - 20Δm,m + 20Δm)
        plt = scatter(meff[range],yerr=Δmeff[range];ylims,label="effective mass")
        add_fit_range!(plt,tmin,tmax,m,Δm;label="")
        plot!(plt,title=title)
        display(plt)
    end
end