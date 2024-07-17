using HiRepParsing
using MixedRepSinglets
using LaTeXStrings
using Plots
include("utils.jl")     
#pgfplotsx(legend=:topright, frame=:box, legendfontsize=14, tickfontsize=14, labelfontsize=14, titlefontsize=16,  markersize=5)
plotlyjs()

path = "/home/fabian/Documents/Physics/Data/DataDiaL/measurementsAS/Lt56Ls24beta6.8mas-1.035/out"
path = "/home/fabian/Downloads/smearing_old_complex_numbers/"
fileCONN = joinpath(path,"out_spectrum_smeared")
fileDISC = joinpath(path,"out_spectrum_smeared_discon")
name  = "F1"
nhits = 128
h5file = "test_F.hdf5"

Nsmear = ("0","40","80")
typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]

typesCONN = r"source_N[0-9]+_sink_N[0-9]+ TRIPLET"
typesDISC = r"DISCON_SEMWALL smear_N[0-9]+ SINGLET"

write = false
if write
    writehdf5_spectrum_with_regexp(fileCONN,h5file,typesCONN,h5group="$name/CONN")
    writehdf5_spectrum_disconnected_with_regexp(fileDISC,h5file,typesDISC,nhits,h5group="$name/DISC")
end

assemble = false
if assemble
    conn, disc = _assemble_correlation_matrix_rep(h5file,name,Nsmear,"";channel="g5",disc_sign=+1,subtract_vev=false,Nf=3)
    correlation_matrix = @. conn - disc
    correlation_matrix_deriv = correlator_derivative(correlation_matrix,t_dim=4)
end

plt3 = plot(legend=:bottomright,legendfontsize=10)
for ind in [(1,1),(3,1),(3,3)] 
    corr = correlation_matrix_deriv[ind[1],ind[2],:,:]
    meff, Δmeff = implicit_meff_jackknife(corr';sign=-1)
    scatter!(plt3,meff, yerr = Δmeff, label="singlet (N1=$(Nsmear[ind[1]]), N2=$(Nsmear[ind[2]]))")
end
eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff_mixed_rep(correlation_matrix;t0=1,binsize=1,deriv=false)
scatter!(plt3,meff[3,:], yerr = Δmeff[3,:],label="singlet (GEVP)")
eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff_mixed_rep(conn;t0=1,binsize=1,deriv=false)
scatter!(plt3,meff[3,:], yerr = Δmeff[3,:],label="non-singlet (GEVP)")
plot!(plt3,ylims=(0,1),xlims=(0,12.5))
display(plt3)