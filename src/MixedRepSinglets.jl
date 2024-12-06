module MixedRepSinglets

    using ProgressMeter
    using Statistics
    using NaNStatistics
    using LinearAlgebra
    using HDF5
    using LatticeUtils
    
    include("smearing_tools.jl")
    export _get_connected_at_smearing_level, _get_disconnected_at_smearing_level
    export _assemble_correlation_matrix_mixed, _assemble_correlation_matrix_rep
    export _assemble_correlation_matrix_rep_nonsinglet
    export stdmean
    
    # Rexports from LatticeUtils
    export eigenvalues, eigenvalues_jackknife_samples, eigenvalues_eigenvectors, eigenvalues_eigenvectors_jackknife_samples
    export implicit_meff
    export add_mass_band!, add_fit_range!, plot_correlator!
    export correlator_derivative
    export awi_corr, awi_fit
    export correlator_folding
    export disconnected_loop_product

end # module MixedRepSinglets
