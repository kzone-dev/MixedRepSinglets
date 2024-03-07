module MixedRepSinglets

    using Statistics
    using NaNStatistics
    using LinearAlgebra
    using Plots
    using LsqFit
    using Roots
    using HDF5

    include("unbiased_estimator.jl")
    export unbiased_estimator, read_hdf5_diagrams
    include("variational_analysis.jl")
    export correlation_matrix, _bin_correlator_matrix, eigenvalues, eigenvalues_jackknife_samples, eigenvalues_eigenvectors
    include("effective_mass.jl")
    export implicit_meff_jackknife, implicit_meff, meff_from_jackknife
    include("plotting.jl")
    export add_mass_band!, plot_correlator!
    include("correlator_derivative.jl")
    export correlator_derivative
    include("fitcorr.jl")
    export fit_corr, fit_corr_bars
    include("pcac.jl")
    export awi_corr, awi_fit
    include("folding.jl")
    export correlator_folding

end # module MixedRepSinglets
