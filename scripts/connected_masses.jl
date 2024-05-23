using DelimitedFiles
using HiRepParsing
using MixedRepSinglets
using HDF5
using Plots
gr(frame=:box)

function main_write_hdf5_logs(path,h5file,parameterfile)
    input = readdlm(parameterfile,';',skipstart=1)
    for prm in eachrow(input)
        
        dir, file, typeCONN, rep, name = prm
        fileCONN = joinpath(path,dir,"out",file)

        @show fileCONN   
        writehdf5_spectrum(fileCONN,h5file,typeCONN,h5group="$name/$rep/CONN")
    end    
end

h5file   = "/home/fabian/Downloads/b55_tests.hdf5"
datapath = "/home/fabian/Dokumente/DataDiaL/"
datapath = "/home/fabian/Documents/DataDiaL/"

isfile(h5file) && rm(h5file)
ispath("output") || mkpath("output")

path = joinpath(datapath,"measurementsTests")
parameterfile655 = "input/parameters_b6p55.csv"
main_write_hdf5_logs(path,h5file,parameterfile655)

path = joinpath(datapath,"measurements")
parameterfile645 = "input/parameters_b6p45.csv"
main_write_hdf5_logs(path,h5file,parameterfile645)

function run_corrfitter(parameters,hdf5file,outdir)
    args = `$(abspath(parameters)) $(abspath(hdf5file)) $(abspath(outdir))`
    try
        run(`python3 $(abspath("scripts/fitting_eigenvalues.py")) $args`)
    catch
        run(`python  $(abspath("scripts/fitting_eigenvalues.py")) $args`)
    end
end

parameters_fitting = "input/parameters_fitting.csv"
run_corrfitter(parameters_fitting,h5file,"output")

results  = readdlm("output/corrfitter_results.csv",';')
fitparam = readdlm("input/parameters_fitting.csv",';',skipstart=1)
for (i,line) in enumerate(eachrow(results))
    ens, channel, rep, T, L, β, m, Δm, χ2dof = line
    tmin, tmax, tp, Nmax = fitparam[i,4:7]
    
    label = "$ens/$rep/CONN/DEFAULT_SEMWALL TRIPLET/$channel"
    corr  = h5read(h5file,label)
    meff, Δmeff = implicit_meff_jackknife(corr')
    #meff = correlator_folding(meff;t_dim=1,sign=+1)

    range = 1:div(T,2)
    ylims = (m - 10Δm,m + 10Δm)
    plt = scatter(meff[range],yerr=Δmeff[range];ylims,label="effective mass")
    display(plt)
end