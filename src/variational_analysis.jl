# choose two sets of operators which shoud lead to identical results. 
# CASE A: η_f: (uΓu + dΓd)/√2 
#         η_a: (UΓU + DΓD + SΓS)/√3
function correlation_matrix(conn_f, conn_as, disc_ff, disc_aa, disc_fa;Nf_fun=2,Nf_as=3)
    #@assert size(conn_f) == size(conn_as) == size(disc_aa) == size(disc_ff) == size(disc_fa) 
    corr_matrix = zeros(eltype(conn_f),(2,2,size(conn_f)...))
    # create the correlators and cross-correlators of CASE 1
    @. corr_matrix[1,1,:,:] = conn_f  - Nf_fun*disc_ff
    @. corr_matrix[2,2,:,:] = conn_as - Nf_as *disc_aa
    @. corr_matrix[1,2,:,:] = sqrt(Nf_fun*Nf_as)*disc_fa
    @. corr_matrix[2,1,:,:] = sqrt(Nf_fun*Nf_as)*disc_fa
    return corr_matrix 
end
function _swap_at_crossing(val,swap)
    N_eig, T = size(val)
    c = swap:T-swap+1
    @assert N_eig == 2
    val[1,c], val[2,c] = val[2,c], val[1,c] 
    return val
end
function eigenvalues(corr;swap=nothing,t0=1)
    eigvals_jk = eigenvalues_jackknife_samples(corr;t0)
    eigvals, Δeigvals = apply_jackknife(eigvals_jk;dims=2)
    if isnothing(swap)
        return eigvals, Δeigvals
    else
        _swap_at_crossing(eigvals, swap)
        _swap_at_crossing(Δeigvals,swap)
        return eigvals, Δeigvals
    end
end
function eigenvalues_eigenvectors(corr;swap=nothing,t0=1)
    eigvals_jk, eigvecs_jk = eigenvalues_eigenvectors_jackknife_samples(corr;t0)
    eigvals, Δeigvals = apply_jackknife(eigvals_jk;dims=2)
    eigvecs, Δeigvecs = apply_jackknife(eigvecs_jk;dims=3)
    if isnothing(swap)
        return eigvals, Δeigvals, eigvecs, Δeigvecs
    else
        _swap_at_crossing(eigvals, swap)
        _swap_at_crossing(Δeigvals,swap)
        # TODO: Sort eigenvectors accordingly
        # TODO: Make it work in this notation
        for i in 1:2
            @views _swap_at_crossing(eigvecs[i,:,:] ,swap)
            @views _swap_at_crossing(Δeigvecs[i,:,:],swap)    
        end
        # Although: The matrix of eigenvectors should be time-independent
        return eigvals, Δeigvals, eigvecs, Δeigvecs
    end
end
function eigenvalues_jackknife_samples(corr;t0 = 1)
    sample = delete1_resample(corr)
    nops, nconf, T = size(sample)[2:4]
    eigvals_jk = zeros(eltype(sample),(nops,nconf,T))
    for s in 1:nconf, t in 1:T
        # smaller values correspond to a faster decay, and thus correspond to a larger masses
        # use sortby to sort the eigenvalues by ascending eigen-energy of the meson state
        eigvals_jk[:,s,t] = eigen(sample[:,:,s,t],sample[:,:,s,t0]).values
    end
    return eigvals_jk
end
function eigenvalues_eigenvectors_jackknife_samples(corr;t0 = 1)
    sample = delete1_resample(corr)
    nops, nconf, T = size(sample)[2:4]
    eigvals_jk = zeros(eltype(sample),(nops,nconf,T))
    eigvecs_jk = zeros(eltype(sample),(nops,nops,nconf,T))
    for s in 1:nconf, t in 1:T
        # smaller values correspond to a faster decay, and thus correspond to a larger masses
        # use sortby to sort the eigenvalues by ascending eigen-energy of the meson state
        sol = eigen(sample[:,:,s,t],sample[:,:,s,t0],sortby= x-> abs(x))
        eigvals_jk[:,s,t] = sol.values
        # I am unsure if the average over all eigenvectors is correct.
        for i in 1:2
            eigvecs_jk[:,i,s,t] = normalize(sol.vectors[:,i])
        end
    end
    return eigvals_jk, eigvecs_jk
end
# generate a resample of the original correlator matrix
function delete1_resample(corr_matrix)
    nops,nconf,T = size(corr_matrix)[2:end]
    samples = similar(corr_matrix)
    # temporary array for jackknife sampling
    tmp = zeros(eltype(corr_matrix),(nops,nops,nconf-1,T))
    for index in 1:nconf    
        for i in 1:index-1
            tmp[:,:,i,:] = corr_matrix[:,:,i,:]
        end
        for i in 1+index:nconf
            tmp[:,:,i-1,:] = corr_matrix[:,:,i,:]
        end
        # perform average after deleting one sample
        samples[:,:,index,:] = dropdims(mean(tmp,dims=3),dims=3)
    end
    return samples
end
# apply jackknife resampling along dimension dims
function apply_jackknife(obs::AbstractArray;dims::Integer)
    N  = size(obs)[dims]
    O  = dropdims(mean(obs;dims);dims)
    ΔO = dropdims(sqrt(N-1)*std(obs;dims,corrected=false);dims)
    return O, ΔO
end
function apply_jackknife(obs::AbstractVector)
    N  = length(obs)
    O  = mean(obs)
    ΔO = sqrt(N-1)*std(obs,corrected=false)
    return O, ΔO
end
