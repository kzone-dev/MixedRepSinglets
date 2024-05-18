function main_write_hdf5_logs(Nsmear,path,hdf5path,parameterfile;filter_channels=true)
    h5file = joinpath(hdf5path,"singlets_smeared.hdf5")

    input = readdlm(parameterfile,';',skipstart=1)
    for prm in eachrow(input)

        dir, patternDISC, patternCONN, nhits, Nmax, rep, name = prm

        fileCONN = joinpath(path,dir,"out/out_spectrum_smeared")
        fileDISC = joinpath(path,dir,"out/out_spectrum_smeared_discon")

        typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
        typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]
        
        channels=["g5","g0g5","g1","g2","g3"]

        @show fileCONN   
        writehdf5_spectrum(fileCONN,h5file,typesCONN,h5group="$name/$rep/CONN";filter_channels,channels)
        @show fileDISC   
        writehdf5_spectrum_disconnected(fileDISC,h5file,typesDISC,nhits,h5group="$name/$rep/DISC";filter_channels,channels)
    end    
end