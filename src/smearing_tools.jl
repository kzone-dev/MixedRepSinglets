function _get_connected_at_smearing_level(h5file,Nsource,Nsink,channel,rep;ensemble="")
    group = "source_N$(Nsource)_sink_N$(Nsink) TRIPLET"
    return h5read(h5file,joinpath(ensemble,rep,"CONN",group,channel))
end
function _get_disconnected_at_smearing_level(h5file,Nsmear,channel,rep;ensemble="")
    group = "DISCON_SEMWALL smear_N$Nsmear SINGLET"
    return h5read(h5file,joinpath(ensemble,rep,"DISC",group,channel))
end
function _assemble_correlation_matrix_mixed(h5file,ensemble,Nsmear;channel="g5",disc_sign=+1,subtract_vev=false)

    discFUN = [_get_disconnected_at_smearing_level(h5file,N,channel,"FUN";ensemble) for N in Nsmear]
    discAS  = [_get_disconnected_at_smearing_level(h5file,N,channel,"AS";ensemble)  for N in Nsmear]
    connFUN = [_get_connected_at_smearing_level(h5file,N1,N2,channel,"FUN";ensemble) for N1 in Nsmear, N2 in Nsmear ]
    connAS  = [_get_connected_at_smearing_level(h5file,N1,N2,channel,"AS";ensemble) for N1 in Nsmear, N2 in Nsmear ]

    # number of configurations and lattice size
    N   = size(first(discFUN))[1]
    T,L = h5read(h5file,joinpath(ensemble,"FUN","CONN","lattice"))[1:2]

    # number of operators in correlation matrix
    S = length(Nsmear)
    Nops = 2*S

    # make sure that there are no NaNs in the data
    @assert 0 == sum(any.(isnan, connFUN))
    @assert 0 == sum(any.(isnan, connAS))
    @assert 0 == sum(any.(isnan, discFUN))
    @assert 0 == sum(any.(isnan, discAS))
    @assert size.(connFUN) == size.(connAS)
    @assert size.(discFUN) == size.(discAS)

    # Compared to the old code, there is another factor of 2 per loop missing
    rescale_disc = 4*L^3
    # rescale connected pieces
    MixedRepSinglets.rescale_connected!.(connFUN ,L)
    MixedRepSinglets.rescale_connected!.(connAS  ,L)

    # model specific parameters
    Nf_fun = 2
    Nf_as  = 3
    
    # create block matrices of the full correlation matrix
    block_diag_FUN = zeros((Nops÷2,Nops÷2,N,T))
    block_diag_AS  = zeros((Nops÷2,Nops÷2,N,T))
    block_mixed    = zeros((Nops÷2,Nops÷2,N,T))
    
    p = Progress( (S^2 + S) ÷ 2 )
    # assemble block matrices for disconnected pieces
    # ( use that the two loops in the disconnected diagra can be interchanged to save computing time) 
    for i in eachindex(Nsmear)
        for j in 1:i
            if i == j
                discFUN_N1N2 = unbiased_estimator(discFUN[i];rescale=rescale_disc,subtract_vev)
                discAS_N1N2  = unbiased_estimator(discAS[i] ;rescale=rescale_disc,subtract_vev)
            else
                discFUN_N1N2 = unbiased_estimator(discFUN[i],discFUN[j];rescale=rescale_disc,subtract_vev) 
                discAS_N1N2  = unbiased_estimator(discAS[i] ,discAS[j] ;rescale=rescale_disc,subtract_vev) 
            end
            discFUNAS_N1N2   = unbiased_estimator(discFUN[i],discAS[j] ;rescale=rescale_disc,subtract_vev) 
            block_diag_FUN[i,j,:,:] = Nf_fun*disc_sign*discFUN_N1N2
            block_diag_FUN[j,i,:,:] = Nf_fun*disc_sign*discFUN_N1N2
            block_diag_AS[i,j,:,:]  = Nf_as *disc_sign*discAS_N1N2
            block_diag_AS[j,i,:,:]  = Nf_as *disc_sign*discAS_N1N2
            block_mixed[i,j,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*discFUNAS_N1N2
            block_mixed[j,i,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*discFUNAS_N1N2
            next!(p) # update progress meter
        end
    end

    # add connected pieces
    for i in eachindex(Nsmear)
        for j in eachindex(Nsmear)
            block_diag_FUN[i,j,:,:] = connFUN[i,j] - block_diag_FUN[i,j,:,:] 
            block_diag_AS[i,j,:,:]  = connAS[i,j]  - block_diag_AS[i,j,:,:]  
        end
    end

    # assemble matrix blocks into full correlation matric
    block_row_1 = vcat(block_diag_FUN,block_mixed)
    block_row_2 = vcat(block_mixed,block_diag_AS)
    correlation_matrix = hcat(block_row_1,block_row_2)

    return correlation_matrix
end
function _assemble_correlation_matrix_rep(h5file,ensemble,Nsmear,rep;channel="g5",disc_sign=+1,subtract_vev=false,Nf)

    disc = [_get_disconnected_at_smearing_level(h5file,N,channel,rep;ensemble) for N in Nsmear]
    conn = [_get_connected_at_smearing_level(h5file,N1,N2,channel,rep;ensemble) for N1 in Nsmear, N2 in Nsmear ]

    # number of configurations and lattice size
    N   = size(first(disc))[1]
    T,L = h5read(h5file,joinpath(ensemble,rep,"CONN","lattice"))[1:2]

    # make sure that there are no NaNs in the data
    @assert 0 == sum(any.(isnan, conn))
    @assert 0 == sum(any.(isnan, disc))

    # Compared to the old code, there is another factor of 2 per loop missing
    rescale_disc = 4*L^3
    # rescale connected pieces
    MixedRepSinglets.rescale_connected!.(conn ,L)

    # number of operators in correlation matrix
    Nops = length(Nsmear)

    # create block matrices of the full correlation matrix
    correlation_matrix_CONN = zeros((Nops,Nops,N,T))
    correlation_matrix_DISC = zeros((Nops,Nops,N,T))
    correlation_matrix_FULL = zeros((Nops,Nops,N,T))
    
    p = Progress( (Nops^2 + Nops) ÷ 2 )
    # assemble block matrices for disconnected pieces
    # ( use that the two loops in the disconnected diagra can be interchanged to save computing time) 
    for i in eachindex(Nsmear)
        for j in 1:i
            if i == j
                disc_N1N2 = unbiased_estimator(disc[i];rescale=rescale_disc,subtract_vev)
            else
                disc_N1N2 = unbiased_estimator(disc[i],disc[j];rescale=rescale_disc,subtract_vev) 
            end
            correlation_matrix_DISC[i,j,:,:] = Nf*disc_sign*disc_N1N2
            correlation_matrix_DISC[j,i,:,:] = Nf*disc_sign*disc_N1N2
            next!(p) # update progress meter
        end
    end

    # add connected pieces
    for i in eachindex(Nsmear)
        for j in eachindex(Nsmear)
            correlation_matrix_CONN[i,j,:,:] = conn[i,j] 
            correlation_matrix_FULL[i,j,:,:] = conn[i,j] - correlation_matrix_DISC[i,j,:,:] 
        end
    end

    return correlation_matrix_CONN, correlation_matrix_DISC
end
function _assemble_correlation_matrix_rep_nonsinglet(h5file,ensemble,Nsmear,rep;channel="g5")

    conn = [_get_connected_at_smearing_level(h5file,N1,N2,channel,rep;ensemble) for N1 in Nsmear, N2 in Nsmear ]

    # number of configurations and lattice size
    N   = size(first(conn))[1]
    T,L = h5read(h5file,joinpath(ensemble,rep,"CONN","lattice"))[1:2]
    MixedRepSinglets.rescale_connected!.(conn ,L)
    Nops = length(Nsmear)

    # make sure that there are no NaNs in the data
    @assert 0 == sum(any.(isnan, conn))
    correlation_matrix_CONN = zeros((Nops,Nops,N,T))
    
    # add connected pieces
    for i in eachindex(Nsmear)
        for j in eachindex(Nsmear)
            correlation_matrix_CONN[i,j,:,:] = conn[i,j] 
        end
    end

    return correlation_matrix_CONN
end
function stdmean(X;dims)
    N = size(X)[dims]
    m = dropdims(mean(X;dims);dims)
    s = dropdims(std(X;dims);dims)/sqrt(N)
    return m, s
end
