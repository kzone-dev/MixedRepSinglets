function run_corrfitter(parameters,hdf5file,outdir)
    args = `$(abspath(parameters)) $(abspath(hdf5file)) $(abspath(outdir))`
    try
        run(`python3 $(abspath("scripts/fitting_eigenvalues.py")) $args`)
    catch
        run(`python  $(abspath("scripts/fitting_eigenvalues.py")) $args`)
    end
end
function main_write_hdf5_logs(path,h5file,parameterfile;regexp=false)
    input = readdlm(parameterfile,';',skipstart=1)
    for prm in eachrow(input)
        
        dir, file, typeCONN, rep, name = prm
        fileCONN = joinpath(path,dir,"out",file)
        @show fileCONN   

        if regexp 
            typeCONN = Regex(typeCONN)
            writehdf5_spectrum_with_regexp(fileCONN,h5file,typeCONN,h5group="$name/$rep/CONN",sort=false)
        else
            writehdf5_spectrum(fileCONN,h5file,typeCONN,h5group="$name/$rep/CONN",sort=false)
        end
        
    end    
end