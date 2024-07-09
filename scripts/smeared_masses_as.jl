using Pkg; Pkg.activate(".")
using DelimitedFiles
using HiRepParsing
using MixedRepSinglets
using HDF5
using Plots
gr(frame=:box)
plotlyjs(frame=:box)
include("utils.jl")

h5file   = "tests_smearing_AS.hdf5"
prm      = "input/parameters_AS.csv"
prm_fit  = "input/fitting_AS.csv"

datapath = "/home/fabian/Dokumente/DataDiaL/"
datapath = "/home/fabian/Documents/Physics/Data/DataDiaL/"

ispath("output") || mkpath("output")
path = joinpath(datapath,"measurementsAS")
Nsmear = 0:30:60

write_hdf5_file = false
if write_hdf5_file
    isfile(h5file) && rm(h5file)
    main_write_hdf5_logs(path,h5file,prm;regexp=true)
end

fitparam = readdlm(prm_fit,';',skipstart=1)
for (i,line) in enumerate(eachrow(fitparam))

    ens, rep, type, channel, tmin, tmax, tp, Nexp  = line
    corr = _assemble_correlation_matrix_rep_nonsinglet(h5file,ens,Nsmear,rep;channel)
    eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff(corr,t0=1)

    
    Nops, T = size(eigvals)

    β   = h5read(h5file,"$ens/$rep/CONN/beta")
    T,L = h5read(h5file,"$ens/$rep/CONN/lattice")[1:2]    
    title = "$T × $L^3, β=$β, $channel, $rep"

    range = 1:div(T,2)
    plt = plot()
    scatter!(plt,meff[Nops,range],yerr=Δmeff[Nops,range];label="effective mass")
    
    @show size(corr)
    for i in eachindex(Nsmear)
        corr0 = corr[i,i,:,:]
        meffd, Δmeffd = implicit_meff_jackknife(corr0')
        scatter!(plt,meffd,yerr=Δmeffd;label="effective mass: diag #$i")
    end

    display(plt)

end
