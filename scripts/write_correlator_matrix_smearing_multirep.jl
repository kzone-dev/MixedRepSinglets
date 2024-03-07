using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using Plots
using HDF5
plotlyjs()
include("smearing_tools.jl")

Nsmear = collect(0:10:80)
h5file = "/home/fabian/Downloads/smeared_singlets_M34.hdf5"
h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M34.hdf5"

for ensemble in ["M3","M4"]

    correlation_matrix = _assemble_correlation_matrix_mixed(h5file,ensemble,Nsmear;channel="g5",disc_sign=+1,subtract_vev=false)

    function _copy_lattice_parameters(outfile,infile,ensemble)
        file = h5open(infile)[joinpath(ensemble,"FUN","CONN")]
        entries = filter(!contains("TRIPLET"),keys(file))
        for entry in entries
            h5write(outfile,joinpath(ensemble,entry),read(file,entry))
        end
    end
    _copy_lattice_parameters(h5corrs,h5file,ensemble)
    h5write(h5corrs,joinpath(ensemble,"correlation_matrix_g5_singlet"),correlation_matrix)
    h5write(h5corrs,joinpath(ensemble,"Wuppertal_levels"),Nsmear)

end