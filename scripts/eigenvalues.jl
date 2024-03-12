function write_eigenvalues(gevp_parameterfile,h5corrs,h5eigenvals)
    parameters = readdlm(gevp_parameterfile,';';skipstart=1)
    for row in eachrow(parameters)

        ensemble, channel, t0, binsize, deriv, ops = row
        nops = parse.(Int,split(replace(ops,r"[()]"=>""),','))
        
        matrixname ="correlation_matrix_$channel"
        correlation_matrix = h5read(h5corrs,joinpath(ensemble,matrixname))
        correlation_matrix = correlation_matrix[nops,nops,:,:]
        
        write_eigenvalues_and_effective_masses(correlation_matrix,h5eigenvals,h5corrs,ensemble,channel;t0,binsize,deriv,resamples=true)
    end
end




