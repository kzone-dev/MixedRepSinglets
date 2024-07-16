using HiRepParsing
using MixedRepSinglets

path = "/home/fabian/Documents/Physics/Data/DataDiaL/measurementsAS"
fileCONN = joinpath(path,"Lt56Ls24beta6.8mas-1.035/out/out_spectrum_smeared")
fileDISC = joinpath(path,"Lt56Ls24beta6.8mas-1.035/out/out_spectrum_smeared_discon")
name  = "AS1"
nhits = 200
h5file = "test_AS.hdf5"


Nsmear = ("0","30","60")
typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]

typesCONN = r"source_N[0-9]+_sink_N[0-9]+ TRIPLET"
typesDISC = r"DISCON_SEMWALL smear_N[0-9]+ SINGLET"

write = false
if write
    writehdf5_spectrum_with_regexp(fileCONN,h5file,typesCONN,h5group="$name/CONN")
    writehdf5_spectrum_disconnected_with_regexp(fileDISC,h5file,typesDISC,nhits,h5group="$name/DISC")
end

re
conn, disc = _assemble_correlation_matrix_rep(h5file,name,Nsmear,"";channel="g5",disc_sign=+1,subtract_vev=false,Nf=3)
correlation_matrix = @. conn - disc

using LaTeXStrings
include("utils.jl")
eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff_mixed_rep(correlation_matrix;t0=1,binsize=1,deriv=true)

using Plots
plotlyjs()
plt = plot()
for i in 3:3
    scatter!(plt,meff[i,:], yerr = Δmeff[i,:],ylims=(0,1),xlims=(0,16))
end
plt