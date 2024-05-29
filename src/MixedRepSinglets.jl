module MixedRepSinglets

    using ProgressMeter
    using Statistics
    using NaNStatistics
    using LinearAlgebra
    using Plots
    using LsqFit
    using Roots
    using HDF5
    using Polyester

    include("unbiased_estimator.jl")
    export unbiased_estimator, read_hdf5_diagrams, unbiased_estimator_threaded
    include("variational_analysis.jl")
    export correlation_matrix, _bin_correlator_matrix, eigenvalues, eigenvalues_jackknife_samples, eigenvalues_eigenvectors, eigenvalues_eigenvectors_jackknife_samples
    include("effective_mass.jl")
    export implicit_meff_jackknife, implicit_meff, meff_from_jackknife
    include("plotting.jl")
    export add_mass_band!, add_fit_range!, plot_correlator!
    include("correlator_derivative.jl")
    export correlator_derivative
    include("fitcorr.jl")
    export fit_corr, fit_corr_bars
    include("pcac.jl")
    export awi_corr, awi_fit
    include("folding.jl")
    export correlator_folding
    include("smearing_tools.jl")
    export _get_connected_at_smearing_level, _get_disconnected_at_smearing_level
    export _assemble_correlation_matrix_mixed, _assemble_correlation_matrix_rep
    export _assemble_correlation_matrix_rep_nonsinglet
    export stdmean

end # module MixedRepSinglets
