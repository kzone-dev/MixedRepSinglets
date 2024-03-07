using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
plotlyjs()

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M34.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M34.hdf5"

function _copy_lattice_parameters(outfile,infile,ensemble)
    file = h5open(infile)[ensemble]
    entries = filter(!contains("correlation_matrix_g5_singlet"),keys(file))
    for entry in entries
        h5write(outfile,joinpath(ensemble,entry),read(file,entry))
    end
end
function eigenvalues_meff_mixed_rep(h5corrs,ensemble;t0 = 1, binsize = 1, deriv = true)

    correlation_matrix = h5read(h5corrs,joinpath(ensemble,"correlation_matrix_g5_singlet"))
    Nsmear = h5read(h5corrs,joinpath(ensemble,"Wuppertal_levels"))
    Nops = length(Nsmear)*2

    symmetry = +1 
    correlation_matrix = correlator_folding(correlation_matrix;t_dim=4,sign=symmetry)
    correlation_matrix = _bin_correlator_matrix(correlation_matrix;binsize)
    if deriv 
        correlation_matrix = correlator_derivative(correlation_matrix;t_dim=4)
        symmetry = -1 
    end

    # use correlator binning
    eigvals, Δeigvals = eigenvalues(correlation_matrix;t0)
    eigenvalues_jackknife = eigenvalues_jackknife_samples(correlation_matrix;t0 ,imag_thresh = 2E-14)
    meff, Δmeff =  meff_from_jackknife(eigenvalues_jackknife;sign=symmetry,swap=nothing)
    return eigvals, Δeigvals, meff, Δmeff
end

for ensemble in ["M3","M4"]

    kws = (t0 = 1, binsize = 1, deriv = true)
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(h5corrs,ensemble;kws...)
    
    write = false
    if write
        _copy_lattice_parameters(h5eigenvals,h5corrs,ensemble)
        h5write(h5eigenvals,joinpath(ensemble,"eigvals"),eigvals)
        h5write(h5eigenvals,joinpath(ensemble,"Delta_eigvals"),Δeigvals)
        h5write(h5eigenvals,joinpath(ensemble,"meff"),meff)
        h5write(h5eigenvals,joinpath(ensemble,"Delta_meff"),Δmeff)
        h5write(h5eigenvals,joinpath(ensemble,"t0"),kws.t0)
        h5write(h5eigenvals,joinpath(ensemble,"deriv"),kws.deriv)
        h5write(h5eigenvals,joinpath(ensemble,"binsize"),kws.binsize)
    end

    range1 = 2:8
    range0 = 2:11
    Nops   = 18
    
    plt = plot()
    scatter!(plt,range0, meff[Nops-0,range0], yerr= Δmeff[Nops-0,range0])
    scatter!(plt,range1, meff[Nops-1,range1], yerr= Δmeff[Nops-1,range1])
    plot!(plt, ylims=(0.25,1))
    display(plt)

end




