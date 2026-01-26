using Pkg; Pkg.activate("."); Pkg.resolve(); Pkg.instantiate(); Pkg.update(); Pkg.precompile() 
Pkg.add("ArgParse")
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--hdf5_folder"
            help = "Path to the directory with the HDF5 files."
            required = true
    end
    return parse_args(s)
end

write_correlator   = true
paramter_path  = "input/parameters/"

args = parse_commandline()
hdf5file_path = args["hdf5_folder"] * "/"

# In order to repsect the dataset size limit on zenodo, only_singlet
# the relevant channels (γ5, γ0γ5, γi) are written to hdf5 file sizes. 
# In order to write all channels to the hdf5 file, set the following
# variable 'write_all_channes_to_hdf5' to 'true'
write_all_channes_to_hdf5 = false

include("run_analysis.jl")  