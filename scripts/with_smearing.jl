using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
plotlyjs()
include("smearing_tools.jl")

Nsmear = 0:10:80
h5file = "/home/fabian/Downloads/smeared_singlets.hdf5"
h5file = "/home/fabian/Downloads/smeared_singlets_M4.hdf5"
io = h5open(h5file)
ensembles = keys(io)
ensemble  = ensembles[1]
close(io)

discFUN = [_get_disconnected_at_smearing_level(h5file,N,"g5","FUN";ensemble) for N in Nsmear]
discAS  = [_get_disconnected_at_smearing_level(h5file,N,"g5","AS";ensemble)  for N in Nsmear]
# first  index: source smearing 
# second index: sink   smearing
connFUN = [_get_connected_at_smearing_level(h5file,N1,N2,"g5","FUN";ensemble) for N1 in Nsmear, N2 in Nsmear ]
connAS  = [_get_connected_at_smearing_level(h5file,N1,N2,"g5","AS";ensemble) for N1 in Nsmear, N2 in Nsmear ]
# choose the smallest value of N for all measurements
N   = size(first(discFUN))[1]
T,L = h5read(h5file,joinpath(ensemble,"FUN","CONN","lattice"))[1:2]

@assert 0 == sum(any.(isnan, connFUN))
@assert 0 == sum(any.(isnan, connAS))
@assert 0 == sum(any.(isnan, discFUN))
@assert 0 == sum(any.(isnan, discAS))

# Compared to the old code, there is another factor of 2 per loop missing
rescale_disc = 4*L^3
# rescale connected pieces
MixedRepSinglets.rescale_connected!.(connFUN ,L)
MixedRepSinglets.rescale_connected!.(connAS  ,L)
# create correlation matrix 
Nf_fun = 2
Nf_as  = 3
disc_sign = +1
subtract_vev = false
Nops = 2*length(Nsmear)

# create block matrices of the full correlation matrix
block_diag_FUN = zeros((Nops÷2,Nops÷2,N,T))
block_diag_AS  = zeros((Nops÷2,Nops÷2,N,T))
block_mixed    = zeros((Nops÷2,Nops÷2,N,T))

# assemble block matrices
for ind1 in eachindex(Nsmear)
    @show Nsmear[ind1]
    for ind2 in eachindex(Nsmear)
        @show Nsmear[ind2]
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
    end
end

block_row_1 = vcat(block_diag_FUN,block_mixed)
block_row_2 = vcat(block_mixed,block_diag_AS)
correlation_matrix = hcat(block_row_1,block_row_2)

plt = plot()
for i in eachindex(Nsmear)
    corr = correlation_matrix[i,i,:,:]
    c, Δc = stdmean(corr,dims=1)
    m, Δm = implicit_meff_jackknife(corr';sign=+1)

    #scatter(c,yerr=Δc,yscale=:log10)
    scatter!(plt,m, yerr= Δm)
end
plt