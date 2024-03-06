using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using ProgressMeter
plotlyjs()
include("smearing_tools.jl")

function _write_smeared_mixed_correlation_matrix_g5(h5file,ensemble,Nsmear)

    discFUN = [_get_disconnected_at_smearing_level(h5file,N,"g5","FUN";ensemble) for N in Nsmear]
    discAS  = [_get_disconnected_at_smearing_level(h5file,N,"g5","AS";ensemble)  for N in Nsmear]
    # first  index: source smearing 
    # second index: sink   smearing
    connFUN = [_get_connected_at_smearing_level(h5file,N1,N2,"g5","FUN";ensemble) for N1 in Nsmear, N2 in Nsmear ]
    connAS  = [_get_connected_at_smearing_level(h5file,N1,N2,"g5","AS";ensemble) for N1 in Nsmear, N2 in Nsmear ]

    # choose the smallest value of N for all measurements
    N   = size(first(discFUN))[1]
    T,L = h5read(h5file,joinpath(ensemble,"FUN","CONN","lattice"))[1:2]

    # make sure that there are no NaNs in the data
    @assert 0 == sum(any.(isnan, connFUN))
    @assert 0 == sum(any.(isnan, connAS))
    @assert 0 == sum(any.(isnan, discFUN))
    @assert 0 == sum(any.(isnan, discAS))

    # Compared to the old code, there is another factor of 2 per loop missing
    rescale_disc = 4*L^3
    # rescale connected pieces
    MixedRepSinglets.rescale_connected!.(connFUN ,L)
    MixedRepSinglets.rescale_connected!.(connAS  ,L)

    # model specific parameters
    Nf_fun = 2
    Nf_as  = 3
    disc_sign = +1
    subtract_vev = false
    
    # number of operators in correlation matrix
    S = length(Nsmear)
    Nops = 2*S

    # create block matrices of the full correlation matrix
    block_diag_FUN = zeros((Nops÷2,Nops÷2,N,T))
    block_diag_AS  = zeros((Nops÷2,Nops÷2,N,T))
    block_mixed    = zeros((Nops÷2,Nops÷2,N,T))
    block_disc_diag_FUN = zeros((Nops÷2,Nops÷2,N,T))
    block_disc_diag_AS  = zeros((Nops÷2,Nops÷2,N,T))
    
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
            block_disc_diag_FUN[i,j,:,:] = Nf_fun*disc_sign*discFUN_N1N2
            block_disc_diag_FUN[j,i,:,:] = Nf_fun*disc_sign*discFUN_N1N2
            block_disc_diag_AS[i,j,:,:]  = Nf_as *disc_sign*discAS_N1N2
            block_disc_diag_AS[j,i,:,:]  = Nf_as *disc_sign*discAS_N1N2
            block_mixed[i,j,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*discFUNAS_N1N2
            block_mixed[j,i,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*discFUNAS_N1N2
            next!(p) # update progress meter
        end
    end

    # add connected pieces
    for i in eachindex(Nsmear)
        for j in eachindex(Nsmear)
            block_diag_FUN[i,j,:,:] = connFUN[i,j] - block_disc_diag_FUN[i,j,:,:] 
            block_diag_AS[i,j,:,:]  = connAS[i,j]  - block_disc_diag_AS[i,j,:,:]  
        end
    end

    # assemble matrix blocks into full correlation matric
    block_row_1 = vcat(block_diag_FUN,block_mixed)
    block_row_2 = vcat(block_mixed,block_diag_AS)
    correlation_matrix = hcat(block_row_1,block_row_2)

    return correlation_matrix
end

Nsmear = 0:10:80
h5file = "/home/fabian/Downloads/smeared_singlets.hdf5"
h5file = "/home/fabian/Downloads/smeared_singlets_M4.hdf5"

correlation_matrix = _write_smeared_correlation_matrix_g5(h5file,"M4",Nsmear)
correlation_matrix_deriv = correlator_derivative(correlation_matrix;t_dim=4)

plt = plot()
for i in eachindex(Nsmear)
    corr = correlation_matrix_deriv[i,i,:,:]
    m, Δm = implicit_meff_jackknife(corr';sign=-1)
    scatter!(plt,m[1:16], yerr= Δm[1:16])
end
plt
