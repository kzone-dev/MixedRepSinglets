using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using ProgressMeter
plotlyjs()
include("smearing_tools.jl")

Nsmear = 0:10:80
h5file = "/home/fabian/Downloads/smeared_singlets_M34.hdf5"

correlation_matrix = _assemble_correlation_matrix_mixed(h5file,"M2",Nsmear)
correlation_matrix_deriv = correlator_derivative(correlation_matrix;t_dim=4)

plt = plot()
for i in eachindex(Nsmear)
    corr = correlation_matrix_deriv[i,i,:,:]
    m, Δm = implicit_meff_jackknife(corr';sign=-1)
    scatter!(plt,m[1:16], yerr= Δm[1:16])
end
plt
