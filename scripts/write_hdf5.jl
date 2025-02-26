function _copy_lattice_parameters(outfile,infile,ensemble;group="")
    file = h5open(infile)[ensemble]
    entries = filter(!contains("correlation_matrix"),keys(file))
    for entry in entries
        label = joinpath(ensemble,group,entry)
        h5write(outfile,label,read(file,entry))
    end
end
function main_write_hdf5_logs(path,hdf5path,parameterfile;filter_channels=true,channels=["g5","g0g5","g1","g2","g3","id"])
    h5file = joinpath(hdf5path,"singlets_smeared.hdf5")

    input = readdlm(parameterfile,';',skipstart=1)
    for prm in eachrow(input)

        dir, typeDISC, typeCONN, fileDISC, fileCONN, nhits, rep, name = prm

        fileCONN = joinpath(path,dir,fileCONN)
        fileDISC = joinpath(path,dir,fileDISC)

        @show dir
        writehdf5_spectrum_with_regexp(fileCONN,h5file,Regex(typeCONN),h5group="$name/$rep/CONN";filter_channels,channels)
        writehdf5_spectrum_disconnected_with_regexp(fileDISC,h5file,Regex(typeDISC),nhits,h5group="$name/$rep/DISC";filter_channels,channels)
    end    
end
function main_write_correlator_matrices(NsmearFUN,NsmearAS,hdf5path)
    h5logfiles = joinpath(hdf5path,"singlets_smeared.hdf5")
    h5corrs    = joinpath(hdf5path,"singlets_smeared_correlators.hdf5")
    isfile(h5corrs) && rm(h5corrs)
    
    # get names of ensembles from hdf5 file
    fid = h5open(h5logfiles, "r")
    ensembles = keys(fid)
    close(fid)

    for ensemble in ensembles

        correlation_matrix_singlet_g5   = _assemble_correlation_matrix_mixed(h5logfiles,ensemble,NsmearFUN,NsmearAS;channel="g5"  ,disc_sign=+1,subtract_vev=false)
        correlation_matrix_singlet_id   = _assemble_correlation_matrix_mixed(h5logfiles,ensemble,NsmearFUN,NsmearAS;channel="id"  ,disc_sign=-1,subtract_vev=false)
        correlation_matrix_singlet_id_vevsubtract = _assemble_correlation_matrix_mixed(h5logfiles,ensemble,NsmearFUN,NsmearAS;channel="id",disc_sign=-1,subtract_vev=true)
        # add non-singlet correlation matrices
        correlation_matrix_nonsinglet_FUN_g5 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearFUN,"FUN";channel="g5")
        correlation_matrix_nonsinglet_FUN_id = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearFUN,"FUN";channel="id")
        correlation_matrix_nonsinglet_AS_g5 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearAS,"AS";channel="g5")
        correlation_matrix_nonsinglet_AS_id = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearAS,"AS";channel="id")
        
        function _copy_lattice_parameters(outfile,infile,ensemble)
            fileFUN = h5open(infile)[joinpath(ensemble,"FUN","CONN")]
            fileAS  = h5open(infile)[joinpath(ensemble,"AS","CONN")]
            
            # ignore everything but correlator
            entries = filter(!contains("TRIPLET"),keys(fileFUN))
            entries = filter(!contains("quarkmasses"),entries)
            for entry in entries
                h5write(outfile,joinpath(ensemble,entry),read(fileFUN,entry))
            end
            # now special case the fermion masses
            h5write(outfile,joinpath(ensemble,"quarkmasses_fundamental")  ,read(fileFUN,"quarkmasses"))
            h5write(outfile,joinpath(ensemble,"quarkmasses_antisymmetric"),read(fileAS, "quarkmasses"))
        end
        _copy_lattice_parameters(h5corrs,h5logfiles,ensemble)
        # NOTE: Note that the entries of the correlation matrix always need to contain the substring "correlation_matrix" 
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g5_singlet"),correlation_matrix_singlet_g5)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_id_singlet"),correlation_matrix_singlet_id)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_id_mvev_singlet"),correlation_matrix_singlet_id_vevsubtract)
        # Add non-singlet correlators here
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g5_nonsinglet_FUN"),correlation_matrix_nonsinglet_FUN_g5)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_id_nonsinglet_FUN"),correlation_matrix_nonsinglet_FUN_id)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g5_nonsinglet_AS"),correlation_matrix_nonsinglet_AS_g5)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_id_nonsinglet_AS"),correlation_matrix_nonsinglet_AS_id)
        # Smearing parameters
        h5write(h5corrs,joinpath(ensemble,"Wuppertal_levels"),NsmearFUN)
        h5write(h5corrs,joinpath(ensemble,"Wuppertal_levels_FUN"),NsmearFUN)
        h5write(h5corrs,joinpath(ensemble,"Wuppertal_levels_AS"),NsmearAS)
        # write single disconnected loops to hdf5 file for analysis of the glueball mixing
        for channel in ["g5","id"] 
            discFUN = stack(_get_disconnected_at_smearing_level(h5logfiles,N,channel,"FUN";ensemble) for N in NsmearFUN)
            discAS  = stack(_get_disconnected_at_smearing_level(h5logfiles,N,channel,"AS";ensemble)  for N in NsmearAS)
            # match the layout of the correlation matrix
            discFUN = permutedims(discFUN,(4,1,2,3))
            discAS  = permutedims(discAS,(4,1,2,3))
            # perform average over sources and use a consistent normalization
            T,L = h5read(h5logfiles,joinpath(ensemble,"FUN","CONN","lattice"))[1:2]
            rescale = sqrt(4*L^3)
            avg_loopFUN = rescale*dropdims(mean(discFUN,dims=3),dims=3)
            avg_loopAS = rescale*dropdims(mean(discAS,dims=3),dims=3)
            # write to HDF5 file
            h5write(h5corrs,joinpath(ensemble,"singlet_loop_$(channel)_FUN"),avg_loopFUN)
            h5write(h5corrs,joinpath(ensemble,"singlet_loop_$(channel)_AS"),avg_loopAS)    
        end
    end
end
function write_eigenvalues(gevp_parameterfile,hdf5path)
    h5corrs     = joinpath(hdf5path,"singlets_smeared_correlators.hdf5")
    h5eigenvals = joinpath(hdf5path,"singlets_smeared_eigenvalues.hdf5")

    isfile(h5eigenvals) && rm(h5eigenvals)

    parameters = readdlm(gevp_parameterfile,';';skipstart=1)
    for row in eachrow(parameters)

        ensemble, channel, t0, binsize, deriv, ops = row
        nops = parse.(Int,split(replace(ops,r"[()]"=>""),','))
        
        matrixname ="correlation_matrix_$channel"
        correlation_matrix = h5read(h5corrs,joinpath(ensemble,matrixname))
        correlation_matrix = correlation_matrix[nops,nops,:,:]
        
        write_eigenvalues_and_effective_masses(correlation_matrix,h5eigenvals,h5corrs,ensemble,channel;t0,binsize,deriv,resamples=true)
    end
end
function write_eigenvalues_and_effective_masses(correlation_matrix,outputfile,inputfile,ensemble,channel; t0, binsize, deriv, resamples = false)
    eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    eigvals_cov = LatticeUtils.cov_jackknife_eigenvalues(eigenvalues_jackknife)
    
    _copy_lattice_parameters(outputfile,inputfile,ensemble;group=channel)

    meff_log, Δmeff_log =  log_meff_jackknife(eigenvalues_jackknife)
    h5write(outputfile,joinpath(ensemble,channel,"meff_log"),meff_log)
    h5write(outputfile,joinpath(ensemble,channel,"Delta_meff_log"),Δmeff_log)

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