using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
using LaTeXStrings
plotlyjs()
pgfplotsx(legend=:topright, frame=:box, legendfontsize=14, tickfontsize=12, labelfontsize=16, markersize=5)

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

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M1.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M1.hdf5"

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M2.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M2.hdf5"

h5corrs = "/home/fabian/Downloads/smeared_singlet_correlators_M34.hdf5"
h5eigenvals = "/home/fabian/Downloads/smeared_singlet_eigenvalues_M34.hdf5"

for ensemble in ["M4"]

    kws = (t0 = 1, binsize = 2, deriv = true)
    eigvals, Δeigvals, meff, Δmeff = eigenvalues_meff_mixed_rep(h5corrs,ensemble;kws...)
    
    write = false
    if write
        _copy_lattice_parameters(h5eigenvals,h5corrs,ensemble)
        h5write(h5eigenvals,joinpath(ensemble,"eigvals_g5_singlet"),eigvals)
        h5write(h5eigenvals,joinpath(ensemble,"Delta_eigvals_g5_singlet"),Δeigvals)
        h5write(h5eigenvals,joinpath(ensemble,"meff_g5_singlet"),meff)
        h5write(h5eigenvals,joinpath(ensemble,"Delta_meff_g5_singlet"),Δmeff)
        h5write(h5eigenvals,joinpath(ensemble,"t0"),kws.t0)
        h5write(h5eigenvals,joinpath(ensemble,"deriv"),kws.deriv)
        h5write(h5eigenvals,joinpath(ensemble,"binsize"),kws.binsize)
    end

    range1 = 2:8
    range0 = 2:11
    Nops   = 18

    β   = h5read(h5corrs,joinpath(ensemble,"beta"))
    T,L = h5read(h5corrs,joinpath(ensemble,"lattice"))[1:2]
    mf  = h5read(h5corrs,joinpath(ensemble,"quarkmasses"))[1]
    mas = "-1.01" #h5read(h5corrs,joinpath(ensemble,"quarkmasses"))[1]
    title = L" N_t \times N_l^3 =%$(T) \times %$(L)^3, \beta=%$β, m_f=%$mf, m_{as}=%$mas"

    plt = plot(title=title, xlabel="t", ylabel="effective mass")
    scatter!(plt,range0, meff[Nops-0,range0], yerr= Δmeff[Nops-0,range0],label=L"$m _{\rm eff}: a $")
    scatter!(plt,range1, meff[Nops-1,range1], yerr= Δmeff[Nops-1,range1],label=L"$m _{\rm eff}: \eta'$")
    plot!(plt, ylims=(0.2,1))
    display(plt)

    isdir("plots/singlet_meff_smeared") || mkdir("plots/singlet_meff_smeared")
    savefig(plt,"plots/singlet_meff_smeared/$(ensemble).pdf")

end




