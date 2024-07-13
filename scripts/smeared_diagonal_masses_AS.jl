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
prm      = "input/parameters_AS.csv"
prm_fit  = "input/fitting_matrix_AS.csv"
datapath = "/home/fabian/Documents/Physics/Data/DataDiaL/"
datapath = "/home/fabian/Dokumente/DataDiaL/"

ispath("output") || mkpath("output")
path = joinpath(datapath,"measurements")
path = joinpath(datapath,"measurementsAS")

write_hdf5_file = true
if write_hdf5_file
    isfile(h5file) && rm(h5file)
    main_write_hdf5_logs(path,h5file,prm;regexp=true)
end

fitparam = readdlm(prm_fit,';',skipstart=1)
#run_corrfitter(prm_fit,h5file,"output")

plt = plot()
for (i,line) in enumerate(eachrow(fitparam))
    ens, rep, type, channel, tmin, tmax, tp, Nexp  = line
    
    β   = h5read(h5file,"$ens/$rep/CONN/beta")
    T,L = h5read(h5file,"$ens/$rep/CONN/lattice")[1:2]
    
    title = "$T × $L^3, β=$β, $channel, $rep"
    
    if channel == "g1"
        label1 = "$ens/$rep/CONN/$type/g1"
        label2 = "$ens/$rep/CONN/$type/g2"
        label3 = "$ens/$rep/CONN/$type/g3"
        corr = (h5read(h5file,label1) .+ h5read(h5file,label2) .+ h5read(h5file,label3)) ./ 3 
    else
        label = "$ens/$rep/CONN/$type/$channel"
        corr = h5read(h5file,label)
    end

    @show size(corr)
    corr = correlator_folding(corr;t_dim=2,sign=+1)
    meff, Δmeff = implicit_meff_jackknife(corr')

    range = 1:div(T,2)
    scatter!(plt,meff[range],yerr=Δmeff[range];label="effective mass ($type)")
    display(plt)

end