using Pkg; Pkg.activate(".")
using DelimitedFiles
using HiRepParsing
using MixedRepSinglets
using HDF5
using Plots
gr(frame=:box)
plotlyjs(frame=:box)
include("utils.jl")

h5file   = "tests_smearing_FUN.hdf5"
prm      = "input/parameters_FUN.csv"
prm_fit  = "input/fitting_FUN.csv"

datapath = "/home/fabian/Dokumente/DataDiaL/"
datapath = "/home/fabian/Documents/Physics/Data/DataDiaL/"

ispath("output") || mkpath("output")
path = joinpath(datapath,"measurementsFUN")
Nsmear = 0:50:100

write_hdf5_file = false
if write_hdf5_file
    isfile(h5file) && rm(h5file)
    main_write_hdf5_logs(path,h5file,prm;regexp=true)
end

fitparam = readdlm(prm_fit,';',skipstart=1)
for (i,line) in enumerate(eachrow(fitparam))

    ens, rep, type, channel, tmin, tmax, tp, Nexp  = line
    corr = _assemble_correlation_matrix_rep_nonsinglet(h5file,ens,Nsmear,rep;channel)
    eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff(corr,t0=2)

    
    Nops, T = size(eigvals)

    β   = h5read(h5file,"$ens/$rep/CONN/beta")
    T,L = h5read(h5file,"$ens/$rep/CONN/lattice")[1:2]    
    title = "$T × $L^3, β=$β, $channel, $rep"

    range = 1:div(T,2)
    plt = plot()
    scatter!(plt,meff[Nops,range],yerr=Δmeff[Nops,range];label="effective mass")
    
    @show size(corr)
    for ind in [(3,1),(2,1)]
        corr0 = corr[ind[1],ind[2],:,:]
        meffd, Δmeffd = implicit_meff_jackknife(corr0')
        scatter!(plt,meffd,yerr=Δmeffd;label="effective mass: C_ij, ij=$ind")
    end

    display(plt)

end
