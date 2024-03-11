using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5

Nsmear = collect(0:10:80)

h5file = "/home/fabian/Downloads/smeared_singlets.hdf5"
h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators.hdf5"

# get names of ensembles from hdf5 file
fid = h5open(h5file, "r")
ensembles = keys(fid)
close(fid)

for ensemble in ensembles

    correlation_matrix = _assemble_correlation_matrix_mixed(h5file,ensemble,Nsmear;channel="g5",disc_sign=+1,subtract_vev=false)

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
    _copy_lattice_parameters(h5corrs,h5file,ensemble)
    h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g5_singlet"),correlation_matrix)
    h5write(h5corrs,joinpath(ensemble,"Wuppertal_levels"),Nsmear)

end