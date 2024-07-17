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

assemble = false
if assemble
    conn, disc = _assemble_correlation_matrix_rep(h5file,name,Nsmear,"";channel="g5",disc_sign=+1,subtract_vev=false,Nf=3)
end
correlation_matrix = @. conn - disc

using LaTeXStrings
using Plots
plotlyjs()
pgfplotsx()
pgfplotsx(legend=:topright, frame=:box, legendfontsize=14, tickfontsize=14, labelfontsize=14, titlefontsize=16,  markersize=5)

include("utils.jl")     
plt = plot(frame=:box)


i = 3
range = 1:13
eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff_mixed_rep(correlation_matrix;t0=1,binsize=1,deriv=false)
scatter!(plt,meff[i,range], yerr = Δmeff[i,range],ylims=(0,1),marker=:circle,label="singlet (GEVP)")
eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff_mixed_rep(conn;t0=1,binsize=1,deriv=false)
scatter!(plt,meff[i,:], yerr = Δmeff[i,:],ylims=(0,1),xlims=(0,22.5),marker=:pent,label="non-singlet (GEVP)")

plotlyjs()
plt3 = plot()
correlation_matrix_deriv = correlator_derivative(conn - disc,t_dim=4)
for ind in [(1,1),(3,1),(3,3)] 
    corrN60N0 = correlation_matrix_deriv[ind[1],ind[2],:,:]
    meffN60N0, ΔmeffN60N0 = implicit_meff_jackknife(corrN60N0';sign=+1)
    scatter!(plt3,meffN60N0, yerr = ΔmeffN60N0,ylims=(0.2,0.6),xlims=(0,22.5),marker=:pent,label="non-singlet (N1=$(Nsmear[ind[1]]), N2=$(Nsmear[ind[2]]))")
end
eigvals, Δeigvals, meff, Δmeff, eigenvalues_jackknife = eigenvalues_meff_mixed_rep(conn;t0=1,binsize=1,deriv=false)
scatter!(plt3,meff[i,:], yerr = Δmeff[i,:],ylims=(0,1),xlims=(0,22.5),marker=:pent,label="non-singlet (GEVP)")
display(plt3)

using Statistics
plotlyjs()
plt2 = plot()
N = size(disc)[3]
disc_avg = dropdims(mean(disc,dims=3),dims=3)
disc_std = dropdims(std(disc,dims=3),dims=3)./sqrt(N)
scatter!(plt2, disc_avg[1,1,:],yerr=disc_std[1,1,:])
scatter!(plt2, disc_avg[2,2,:],yerr=disc_std[2,2,:])
scatter!(plt2, disc_avg[3,3,:],yerr=disc_std[3,3,:])

plt
plt2
2+2