using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using HDF5

function _fold_corr!(corr;sign=+1)
    nop, nop, N, T = size(corr)
    # Skip first entry since it does not to be averaged
    # Add 1 to all indices: julia is one-indexed
    for t in 1:div(T,2)
        t1 = t+1
        t2 = T-t+1
        tmp1 = (corr[:,:,:,t1] + sign*corr[:,:,:,t2])/2
        tmp2 = (sign*corr[:,:,:,t1] + corr[:,:,:,t2])/2
        corr[:,:,:,t1] = tmp1
        corr[:,:,:,t2] = tmp2
    end
end
function writehdf5_correlationmatrix(h5file,fileCONN_fun,fileDISC_fun,fileCONN_as,fileDISC_as;Γ="g5")
    HiRepParsing._write_lattice_setup(fileCONN_fun,h5file)
    
    conn_f, conn_as, disc_ff, disc_aa, disc_fa = read_hdf5_diagrams(fileCONN_fun,fileDISC_fun,fileCONN_as,fileDISC_as;Γ)
    corr = correlation_matrix(conn_f,conn_as,disc_ff,disc_aa,disc_fa)
    
    h5write(h5file,"singlet_correlation_matrix_$Γ",corr)
    _fold_corr!(corr;sign=+1)
    h5write(h5file,"singlet_correlation_matrix_$(Γ)_folded",corr)
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
ensembles = [
    "Lt48Ls20beta6.5mf0.71mas1.01",
    "Lt64Ls20beta6.5mf0.71mas1.01",
    "Lt64Ls20beta6.5mf0.70mas1.01",
    "Lt80Ls20beta6.5mf0.71mas1.01", 
    "Lt96Ls20beta6.5mf0.71mas1.01"
]

for name in ensembles
    writehdf5_correlationmatrix(name,savedir="output/correlation_matrix")
end