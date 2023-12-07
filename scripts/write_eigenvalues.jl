using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HDF5
using HiRepParsing

# In these correlators the periodicity is such that c[2] = c[T]
basepath = "./output/correlation_matrix/"
savepath = "./output/eigenvalues/"

ensembles = [
    "Lt48Ls20beta6.5mf0.71mas1.01",
    "Lt64Ls20beta6.5mf0.71mas1.01",
    "Lt64Ls20beta6.5mf0.70mas1.01",
    "Lt80Ls20beta6.5mf0.71mas1.01", 
    "Lt96Ls20beta6.5mf0.71mas1.01"
]
swaps = [6, 6, 6, 6, 6]

for (i,name) in enumerate(ensembles)
    Γ    = "g5"
    file = joinpath(basepath,"correlation_matrix_$name.h5")
    h5file = joinpath(savepath,"eigenvalues_$name.h5")

    corr = h5read(file,"singlet_correlation_matrix_$Γ")
    corr = _bin_correlator_matrix(corr;binsize=2) 
    corr_folded = h5read(file,"singlet_correlation_matrix_$(Γ)_folded")
    corr_folded = _bin_correlator_matrix(corr_folded;binsize=2) 

    ispath(savepath) || mkpath(savepath)

    # write lattice parameters
    for key in ["plaquette","configurations","gauge group","quarkmasses","beta","lattice"]
        h5write(h5file,key,h5read(file,key))
    end
   
    # get eigenvalues
    ev, Δev, vecs, Δvecs = eigenvalues_eigenvectors(corr;swap=swaps[i],t0=1)
    evF, ΔevF, vecsF, ΔvecsF = eigenvalues_eigenvectors(corr_folded;swap=swaps[i],t0=1)
    h5write(h5file,"singlet_eigenvalues_$Γ",ev)
    h5write(h5file,"singlet_eigenvectors_$Γ",vecs)
    h5write(h5file,"Delta_singlet_eigenvalues_$Γ",Δev)
    h5write(h5file,"Delta_singlet_eigenvectors_$Γ",Δvecs)
    h5write(h5file,"singlet_eigenvalues_$(Γ)_folded",evF)
    h5write(h5file,"singlet_eigenvectors_$(Γ)_folded",vecsF)
    h5write(h5file,"Delta_singlet_eigenvalues_$(Γ)_folded",ΔevF)
    h5write(h5file,"Delta_singlet_eigenvectors_$(Γ)_folded",ΔvecsF)
end