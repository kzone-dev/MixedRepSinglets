using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5

function writehdf5_correlationmatrix(h5file,fileCONN_fun,fileDISC_fun,fileCONN_as,fileDISC_as;Γ="g5")
    HiRepParsing._write_lattice_setup(fileCONN_fun,h5file)
    
    conn_f, conn_as, disc_ff, disc_aa, disc_fa = read_hdf5_diagrams(fileCONN_fun,fileDISC_fun,fileCONN_as,fileDISC_as;Γ)
    corr = correlation_matrix(conn_f,conn_as,disc_ff,disc_aa,disc_fa)
    
    h5write(h5file,"singlet_correlation_matrix_$Γ",corr)   
end
function writehdf5_correlationmatrix(ensemble_name;savedir="",kws...)
    fileCONN_fun = joinpath(basepath,"$(ensemble_name)FUN/out_spectrum.h5")
    fileCONN_as  = joinpath(basepath,"$(ensemble_name)AS/out_spectrum.h5")
    fileDISC_fun = joinpath(basepath,"$(ensemble_name)FUN/out_spectrum_discon.h5")
    fileDISC_as  = joinpath(basepath,"$(ensemble_name)AS/out_spectrum_discon.h5")

    ispath(savedir) || mkpath(savedir)
    h5file = joinpath(savedir,"correlation_matrix_$ensemble_name.h5")
    writehdf5_correlationmatrix(h5file,fileCONN_fun,fileDISC_fun,fileCONN_as,fileDISC_as;kws...)
end

basepath = "./output/h5files/"
ensembles = ["Lt64Ls20beta6.5mf0.71mas1.01", "Lt64Ls20beta6.5mf0.70mas1.01", "Lt48Ls20beta6.5mf0.71mas1.01"]
for name in ensembles
    writehdf5_correlationmatrix(name,savedir="output/correlation_matrix")
end