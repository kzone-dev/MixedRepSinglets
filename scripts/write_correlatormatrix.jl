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
        correlation_matrix_singlet_g0g5 = _assemble_correlation_matrix_mixed(h5logfiles,ensemble,NsmearFUN,NsmearAS;channel="g0g5",disc_sign=+1,subtract_vev=false)
        correlation_matrix_singlet_id   = _assemble_correlation_matrix_mixed(h5logfiles,ensemble,NsmearFUN,NsmearAS;channel="id"  ,disc_sign=-1,subtract_vev=false)
        correlation_matrix_singlet_id_vevsubtract = _assemble_correlation_matrix_mixed(h5logfiles,ensemble,NsmearFUN,NsmearAS;channel="id",disc_sign=-1,subtract_vev=true)
        correlation_matrix_nonsinglet_FUN_g5 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearFUN,"FUN";channel="g5")
        correlation_matrix_nonsinglet_FUN_g1 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearFUN,"FUN";channel="g1")
        correlation_matrix_nonsinglet_FUN_g2 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearFUN,"FUN";channel="g2")
        correlation_matrix_nonsinglet_FUN_g3 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearFUN,"FUN";channel="g3")
        correlation_matrix_nonsinglet_AS_g5 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearAS,"AS";channel="g5")
        correlation_matrix_nonsinglet_AS_g1 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearAS,"AS";channel="g1")
        correlation_matrix_nonsinglet_AS_g2 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearAS,"AS";channel="g2")
        correlation_matrix_nonsinglet_AS_g3 = _assemble_correlation_matrix_rep_nonsinglet(h5logfiles,ensemble,NsmearAS,"AS";channel="g3")

        correlation_matrix_nonsinglet_FUN_g1 = @. (correlation_matrix_nonsinglet_FUN_g1 + correlation_matrix_nonsinglet_FUN_g2 + correlation_matrix_nonsinglet_FUN_g3)/3
        correlation_matrix_nonsinglet_AS_g1  = @. (correlation_matrix_nonsinglet_AS_g1  + correlation_matrix_nonsinglet_AS_g2  + correlation_matrix_nonsinglet_AS_g3 )/3

        function _copy_lattice_parameters(outfile,infile,ensemble)
            fileFUN = h5open(infile)[joinpath(ensemble,"FUN","CONN")]
            fileAS  = h5open(infile)[joinpath(ensemble,"AS","CONN")]
            
            # ignore everything but correlator
            entries = filter(!contains("TRIPLET"),keys(fileFUN))
            entries = filter(!contains("quarkmasses"),entries)
            @show entries
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
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g0g5_singlet"),correlation_matrix_singlet_g0g5)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_id_singlet"),correlation_matrix_singlet_id)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_id_mvev_singlet"),correlation_matrix_singlet_id_vevsubtract)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g5_nonsinglet_FUN"),correlation_matrix_nonsinglet_FUN_g5)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g1_nonsinglet_FUN"),correlation_matrix_nonsinglet_FUN_g1)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g5_nonsinglet_AS"),correlation_matrix_nonsinglet_AS_g5)
        h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g1_nonsinglet_AS"),correlation_matrix_nonsinglet_AS_g1)
        h5write(h5corrs,joinpath(ensemble,"Wuppertal_levels"),NsmearFUN)
        h5write(h5corrs,joinpath(ensemble,"Wuppertal_levels_FUN"),NsmearFUN)
        h5write(h5corrs,joinpath(ensemble,"Wuppertal_levels_AS"),NsmearAS)
    end
end