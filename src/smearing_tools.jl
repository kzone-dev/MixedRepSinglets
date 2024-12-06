function _get_connected_at_smearing_level(h5file,Nsource,Nsink,channel,rep;ensemble="")
    group = "source_N$(Nsource)_sink_N$(Nsink) TRIPLET"
    data  = h5read(h5file,joinpath(ensemble,rep,"CONN",group,channel))
    @assert 0 == sum(any.(isnan, data))
    return data 
end
function _get_disconnected_at_smearing_level(h5file,Nsmear,channel,rep;ensemble="")
    group = "DISCON_SEMWALL smear_N$Nsmear SINGLET"
    data  =  h5read(h5file,joinpath(ensemble,rep,"DISC",group,channel))
    @assert 0 == sum(any.(isnan, data))
    return data 
end
function rescale_connected!(corr,L)
    n1 = L^6/2 # from the norm used in HiReo
    n2 = L^3   # only keep a norm for the spatial volume
    @. corr  *= (n1/n2)
end
_assemble_correlation_matrix_mixed(h5file,ensemble,Nsmear;kws...) = _assemble_correlation_matrix_mixed(h5file,ensemble,Nsmear,Nsmear;kws...)
function _assemble_correlation_matrix_mixed(h5file,ensemble,NsmearFUN,NsmearAS;channel="g5",disc_sign=+1,subtract_vev=false,nsrc_max=typemax(Int64))

    discFUN = [_get_disconnected_at_smearing_level(h5file,N,channel,"FUN";ensemble) for N in NsmearFUN]
    discAS  = [_get_disconnected_at_smearing_level(h5file,N,channel,"AS";ensemble)  for N in NsmearAS]
    connFUN = [_get_connected_at_smearing_level(h5file,N1,N2,channel,"FUN";ensemble) for N1 in NsmearFUN, N2 in NsmearFUN ]
    connAS  = [_get_connected_at_smearing_level(h5file,N1,N2,channel,"AS";ensemble) for N1 in NsmearAS, N2 in NsmearAS ]

    # number of configurations and lattice size
    N   = size(first(discFUN))[1]
    T,L = h5read(h5file,joinpath(ensemble,"FUN","CONN","lattice"))[1:2]
    NF  = length(NsmearFUN)
    NA  = length(NsmearAS)

    # Compared to the old code, there is another factor of 2 per loop missing
    rescale_disc = 4*L^3
    rescale_connected!.(connFUN ,L) # rescale connected pieces
    rescale_connected!.(connAS  ,L) # rescale connected pieces

    # model specific parameters
    Nf_fun = 2
    Nf_as  = 3

    correlation_matrix = zeros((NF+NA,NF+NA,N,T))
    p = Progress( (NF^2 + NF) ÷ 2 + (NA^2 + NA) ÷ 2 + NF*NA )
    for i in eachindex(NsmearFUN)
        for j in 1:i
            if i == j
                discFUN_N1N2 = disconnected_loop_product(discFUN[i];rescale=rescale_disc,subtract_vev,nsrc_max)
            else
                discFUN_N1N2 = disconnected_loop_product(discFUN[i],discFUN[j];rescale=rescale_disc,subtract_vev,nsrc_max) 
            end
            correlation_matrix[i,j,:,:] = connFUN[i,j] - Nf_fun*disc_sign*discFUN_N1N2
            correlation_matrix[j,i,:,:] = connFUN[j,i] - Nf_fun*disc_sign*discFUN_N1N2
            next!(p) # update progress meter
        end
    end
    for i in eachindex(NsmearAS)
        for j in 1:i
            if i == j
                discAS_N1N2  = disconnected_loop_product(discAS[i] ;rescale=rescale_disc,subtract_vev,nsrc_max)
            else
                discAS_N1N2  = disconnected_loop_product(discAS[i] ,discAS[j] ;rescale=rescale_disc,subtract_vev,nsrc_max) 
            end
            correlation_matrix[i+NF,j+NF,:,:] = connAS[i,j] - Nf_as *disc_sign*discAS_N1N2
            correlation_matrix[j+NF,i+NF,:,:] = connAS[j,i] - Nf_as *disc_sign*discAS_N1N2
            next!(p) # update progress meter
        end
    end
    for i in eachindex(NsmearFUN)
        for j in eachindex(NsmearAS)
            discFUNAS = disconnected_loop_product(discFUN[i],discAS[j] ;rescale=rescale_disc,subtract_vev,nsrc_max) 
            correlation_matrix[i,j+NF,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*discFUNAS
            correlation_matrix[j+NF,i,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*discFUNAS
            next!(p) # update progress meter
        end
    end
    return correlation_matrix
end
function _assemble_correlation_matrix_rep(h5file,ensemble,Nsmear,rep;channel="g5",disc_sign=+1,subtract_vev=false,Nf,nsrc_max=typemax(Int64))

    disc = [_get_disconnected_at_smearing_level(h5file,N,channel,rep;ensemble) for N in Nsmear]
    conn = [_get_connected_at_smearing_level(h5file,N1,N2,channel,rep;ensemble) for N1 in Nsmear, N2 in Nsmear ]

    # number of configurations and lattice size
    N    = size(first(disc))[1]
    T,L  = h5read(h5file,joinpath(ensemble,rep,"CONN","lattice"))[1:2]
    Nops = length(Nsmear)

    # Compared to the old code, there is another factor of 2 per loop missing
    rescale_disc = 4*L^3
    rescale_connected!.(conn ,L) # rescale connected pieces

    # create block matrices of the full correlation matrix
    correlation_matrix_CONN = zeros((Nops,Nops,N,T))
    correlation_matrix_DISC = zeros((Nops,Nops,N,T))
    
    p = Progress( (Nops^2 + Nops) ÷ 2 )
    for i in eachindex(Nsmear)
        for j in 1:i
            if i == j
                disc_N1N2 = disconnected_loop_product(disc[i];rescale=rescale_disc,subtract_vev,nsrc_max)
            else
                disc_N1N2 = disconnected_loop_product(disc[i],disc[j];rescale=rescale_disc,subtract_vev,nsrc_max) 
            end
            correlation_matrix_DISC[i,j,:,:] = Nf*disc_sign*disc_N1N2
            correlation_matrix_DISC[j,i,:,:] = Nf*disc_sign*disc_N1N2
            correlation_matrix_CONN[i,j,:,:] = conn[i,j] 
            correlation_matrix_CONN[j,i,:,:] = conn[j,i] 
            next!(p) # update progress meter
        end
    end
    return correlation_matrix_CONN, correlation_matrix_DISC
end
function _assemble_correlation_matrix_rep_nonsinglet(h5file,ensemble,Nsmear,rep;channel="g5")
    conn = [_get_connected_at_smearing_level(h5file,N1,N2,channel,rep;ensemble) for N1 in Nsmear, N2 in Nsmear ]
    # number of configurations and lattice size
    N   = size(first(conn))[1]
    T,L = h5read(h5file,joinpath(ensemble,rep,"CONN","lattice"))[1:2]
    rescale_connected!.(conn ,L)
    Nops = length(Nsmear)
    correlation_matrix_CONN = zeros((Nops,Nops,N,T))    
    # add connected pieces
    for i in eachindex(Nsmear),j in eachindex(Nsmear)
        correlation_matrix_CONN[i,j,:,:] = conn[i,j] 
    end
    return correlation_matrix_CONN
end
function stdmean(X;dims)
    N = size(X)[dims]
    m = dropdims(mean(X;dims);dims)
    s = dropdims(std(X;dims);dims)/sqrt(N)
    return m, s
end
