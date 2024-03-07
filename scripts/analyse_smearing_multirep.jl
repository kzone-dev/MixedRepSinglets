using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
plotlyjs()

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M34.hdf5"
ensemble = "M3"

correlation_matrix = h5read(h5corrs,joinpath(ensemble,"correlation_matrix_g5_singlet"))
Nsmear = h5read(h5corrs,joinpath(ensemble,"Wuppertal_levels"))
Nops = length(Nsmear)*2

symmetry = +1 
correlation_matrix = correlator_folding(correlation_matrix;t_dim=4,sign=symmetry)
correlation_matrix = _bin_correlator_matrix(correlation_matrix;binsize=2)
correlation_matrix_deriv = correlator_derivative(correlation_matrix;t_dim=4)
symmetry = -1 

t0 = 1
# use correlator binning
eigvals, Δeigvals = eigenvalues(correlation_matrix_deriv;t0)
eigenvalues_jackknife = eigenvalues_jackknife_samples(correlation_matrix_deriv;t0 ,imag_thresh = 2E-14)
meff, Δmeff =  meff_from_jackknife(eigenvalues_jackknife;sign=symmetry,swap=nothing)

range1 = 2:8
range0 = 2:11

plt = plot()
scatter!(plt,range0, meff[Nops-0,range0], yerr= Δmeff[Nops-0,range0])
scatter!(plt,range1, meff[Nops-1,range1], yerr= Δmeff[Nops-1,range1])
plot!(ylims=(0.25,1))
