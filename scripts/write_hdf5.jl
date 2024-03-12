function main_write_hdf5_logs(Nsmear,path,h5file,parameterfile)
    input = readdlm(parameterfile,';',skipstart=1)
    for prm in eachrow(input)

        dir, patternDISC, patternCONN, nhits, Nmax, rep, name = prm

        fileCONN = joinpath(path,dir,"out/out_spectrum_smeared")
        fileDISC = joinpath(path,dir,"out/out_spectrum_smeared_discon")

        typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
        typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]
        
        @show fileCONN   
        writehdf5_spectrum(fileCONN,h5file,typesCONN,h5group="$name/$rep/CONN")
        @show fileDISC   
        writehdf5_spectrum_disconnected(fileDISC,h5file,typesDISC,nhits,h5group="$name/$rep/DISC")
    end    
end