using Pkg; Pkg.activate(".")
using DelimitedFiles
using HiRepParsing
using MixedRepSinglets
using HDF5
using Plots
gr(frame=:box)
plotlyjs(frame=:box)
include("utils.jl")

h5file   = "/home/fabian/Downloads/tests_smearing.hdf5"
prm      = "input/parameters_smearing.csv"
prm_fit  = "input/parameters_fitting_smearing.csv"
datapath = "/home/fabian/Dokumente/DataDiaL/"
datapath = "/home/fabian/Documents/Physics/Data/DataDiaL/"

ispath("output") || mkpath("output")
path = joinpath(datapath,"measurements")
Nsmear = 0:40:120

write_hdf5_file = true
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
    plot!(plt,title=title,ylims=(0.34,0.36))
    display(plt)

end