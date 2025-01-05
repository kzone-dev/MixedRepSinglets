start_from_logs    = false
write_correlator   = false
write_gevp_results = true

logfiles_path  = "/home/fabian/Dokumente/Physics/Data/DataDiaL/measurements/"
logfiles_path  = "/home/fabian/Documents/Physics/Data/DataDiaL/measurements/"
paramter_path  = "input/parameters_M3/"
hdf5file_path  = "output/hdf5out"
output_path    = "output/"

# In order to repsect the dataset size limit on zenodo, only_singlet
# the relevant channels (γ5, γ0γ5, γi) are written to hdf5 file sizes. 
# In order to write all channels to the hdf5 file, set the following
# variable 'write_all_channes_to_hdf5' to 'true'
write_all_channes_to_hdf5 = false

include("run_analysis.jl")  