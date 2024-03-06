using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using ProgressMeter
plotlyjs()
include("smearing_tools.jl")

function _write_smeared_correlation_matrix_g5(h5file,ensemble,Nsmear)

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
    Nops = 2*length(Nsmear)

    # create block matrices of the full correlation matrix
    block_diag_FUN = zeros((Nops÷2,Nops÷2,N,T))
    block_diag_AS  = zeros((Nops÷2,Nops÷2,N,T))
    block_mixed    = zeros((Nops÷2,Nops÷2,N,T))

    p = Progress(length(Nsmear)^2)
    # assemble block matrices
    for ind1 in eachindex(Nsmear)
        for ind2 in eachindex(Nsmear)
            if ind1 == ind2
                discFUN_N1N2 = unbiased_estimator(discFUN[ind1];rescale=rescale_disc,subtract_vev)
                discAS_N1N2  = unbiased_estimator(discAS[ind1] ;rescale=rescale_disc,subtract_vev)
            else
                discFUN_N1N2 = unbiased_estimator(discFUN[ind1],discFUN[ind2];rescale=rescale_disc,subtract_vev) 
                discAS_N1N2  = unbiased_estimator(discAS[ind1] ,discAS[ind2] ;rescale=rescale_disc,subtract_vev) 
            end
            discFUNAS_N1N2   = unbiased_estimator(discFUN[ind1],discAS[ind2] ;rescale=rescale_disc,subtract_vev) 
            block_diag_FUN[ind1,ind2,:,:] = connFUN[ind1,ind2] - Nf_fun*disc_sign*discFUN_N1N2
            block_diag_AS[ind1,ind2,:,:]  = connAS[ind1,ind2]  - Nf_as *disc_sign*discAS_N1N2
            block_mixed[ind1,ind2,:,:] = sqrt(Nf_fun*Nf_as)*disc_sign*discFUNAS_N1N2
            next!(p)
        end
    end

    block_row_1 = vcat(block_diag_FUN,block_mixed)
    block_row_2 = vcat(block_mixed,block_diag_AS)

    correlation_matrix = hcat(block_row_1,block_row_2)
    return correlation_matrix
end

Nsmear = 0:40:80
h5file = "/home/fabian/Downloads/smeared_singlets.hdf5"
h5file = "/home/fabian/Downloads/smeared_singlets_M4.hdf5"

correlation_matrix = _write_smeared_correlation_matrix_g5(h5file,ensemble,Nsmear)
@profview _write_smeared_correlation_matrix_g5(h5file,ensemble,Nsmear)
#correlation_matrix = _write_smeared_correlation_matrix_g5(h5file,ensemble,Nsmear)
#correlation_matrix_deriv = correlator_derivative(correlation_matrix;t_dim=4)

plt = plot()
for i in 1:3 #eachindex(Nsmear)
    corr = correlation_matrix[i,i,:,:]
    corr = correlation_matrix_deriv[i,i,:,:]
    c, Δc = stdmean(corr,dims=1)
    m, Δm = implicit_meff_jackknife(corr';sign=-1)

    #scatter(c,yerr=Δc,yscale=:log10)
    scatter!(plt,m[1:16], yerr= Δm[1:16])
end
plt
