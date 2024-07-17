using HiRepParsing
using MixedRepSinglets
using LaTeXStrings
using Plots
include("utils.jl")     
#pgfplotsx(legend=:topright, frame=:box, legendfontsize=14, tickfontsize=14, labelfontsize=14, titlefontsize=16,  markersize=5)
gr(frame=:box,legendfontsize=8,legend=:bottomright)

path = "/home/fabian/Downloads/smearing_old_complex_numbers/"
fileCONN = joinpath(path,"out_spectrum_smeared")
fileDISC = joinpath(path,"out_spectrum_smeared_discon")
name  = "F1"
nhits = 128
h5file = "test_F.hdf5"
Nsmear = ("0","40","80")
mπREF, ΔmπREF = 0.5663, 0.0008  
mηREF, ΔmηREF = 0.6100, 0.0060  

path = "/home/fabian/Documents/Physics/Data/DataDiaL/measurementsAS/Lt56Ls24beta6.8mas-1.035/out"
fileCONN = joinpath(path,"out_spectrum_smeared")
fileDISC = joinpath(path,"out_spectrum_smeared_discon")
name  = "AS1"
nhits = 200
h5file = "test_AS.hdf5"
Nsmear = ("0","30","60")
mπREF, ΔmπREF = NaN, NaN
mηREF, ΔmηREF = NaN, NaN

typesDISC = ["DISCON_SEMWALL smear_N$N SINGLET"  for N  in Nsmear]
typesCONN = ["source_N$(N1)_sink_N$(N2) TRIPLET" for N1 in Nsmear, N2 in Nsmear]

typesCONN = r"source_N[0-9]+_sink_N[0-9]+ TRIPLET"
typesDISC = r"DISCON_SEMWALL smear_N[0-9]+ SINGLET"

write = false
if write
    writehdf5_spectrum_with_regexp(fileCONN,h5file,typesCONN,h5group="$name/CONN")
    writehdf5_spectrum_disconnected_with_regexp(fileDISC,h5file,typesDISC,nhits,h5group="$name/DISC")
end

assemble = true
if assemble
    conn, disc = _assemble_correlation_matrix_rep(h5file,name,Nsmear,"";channel="g5",disc_sign=+1,subtract_vev=false,Nf=3)
    correlation_matrix = @. conn - disc
    correlation_matrix_deriv = correlator_derivative(correlation_matrix,t_dim=4)
end

binsize = 1
deriv = false
t0 = 1
Nl = length(Nsmear)

plt3 = plot()
for ind in [(1,1),(Nl,1),(Nl,Nl)] 
    corr = deriv ? correlation_matrix_deriv[ind[1],ind[2],:,:] : correlation_matrix[ind[1],ind[2],:,:]
    sign = deriv ? -1 : +1
    meff, Δmeff = implicit_meff_jackknife(corr';sign)
    label="singlet (N1=$(Nsmear[ind[1]]), N2=$(Nsmear[ind[2]]))"
    plot!(plt3,meff, yerr = Δmeff, ms=5, markershape=:auto; label)
end
eigvals1, Δeigvals1, meff1, Δmeff1, eigenvalues_jackknife1 = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
eigvals2, Δeigvals2, meff2, Δmeff2, eigenvalues_jackknife2 = eigenvalues_meff_mixed_rep(conn;t0,binsize,deriv=false)
plot!(plt3,meff1[Nl,:], yerr = Δmeff1[Nl,:],label="singlet (GEVP)",ms=5,markershape=:auto)
plot!(plt3,meff2[Nl,:], yerr = Δmeff2[Nl,:],label="non-singlet (GEVP)",ms=5,markershape=:auto)
#add_mass_band!(plt3,mπREF,ΔmπREF;label="",alpha=0.5)
#add_mass_band!(plt3,mηREF,ΔmηREF;label="",alpha=0.5)
plot!(plt3,ylims=(0,1),xlims=(0,12.5),legend=:bottomleft)
display(plt3)