function main_write_hdf5_logs(path,hdf5path,parameterfile;filter_channels=true)
    h5file = joinpath(hdf5path,"singlets_smeared.hdf5")

    input = readdlm(parameterfile,';',skipstart=1)
    for prm in eachrow(input)

        dir, typeCONN, typeDISC, fileCONN, fileDISC, nhits, rep, name = prm

        fileCONN = joinpath(path,dir,fileCONN)
        fileDISC = joinpath(path,dir,fileDISC)
        channels=["g5","g0g5","g1","g2","g3","id"]

        @show dir
        writehdf5_spectrum_with_regexp(fileCONN,h5file,Regex(typeCONN),h5group="$name/$rep/CONN";filter_channels,channels)
        writehdf5_spectrum_disconnected_with_regexp(fileDISC,h5file,Regex(typeDISC),nhits,h5group="$name/$rep/DISC";filter_channels,channels)
    end    
end