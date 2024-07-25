using HiRepParsing
using MixedRepSinglets
using LaTeXStrings
using Plots
include("utils.jl")
include("utils_single_rep.jl")
function plot_from_logs(fileCONN,fileDISC,h5file,name,Nsmear,nhits;write,scalar,Nf,xlim,Ns_conn,plot_sing,plot_conn,deriv)
    meff1, Δmeff1, meff2, Δmeff2, conn, correlation_matrix, correlation_matrix_deriv = correlation_matrices_single_rep(fileCONN,fileDISC,h5file,name,Nsmear,nhits;write,scalar,Nf,deriv)
    plt = plot_singlet_vs_nonsinglets(meff1, Δmeff1, meff2, Δmeff2, conn, correlation_matrix, correlation_matrix_deriv, Nsmear; xlim, Ns_conn, deriv, plot_sing, plot_conn)
    return plt
end
gr(frame=:box,legendfontsize=8,legend=:bottomright)
plotlyjs(frame=:box,legendfontsize=8,legend=:bottomright)
pgfplotsx(legendfonthalign = :left,legend=:topright, frame=:box, legendfontsize=12, tickfontsize=14, labelfontsize=14, titlefontsize=16,  markersize=5)

path = "/home/fabian/Downloads/smearing_old_complex_numbers/"
fileCONN = joinpath(path,"out_spectrum_smeared")
fileDISC = joinpath(path,"out_spectrum_smeared_discon")
mπREF, ΔmπREF = 0.5663, 0.0008  
mVREF, ΔmVREF = 0.6994, 0.0018
mηREF, ΔmηREF = 0.6100, 0.0060
h5file = "test.hdf5"
name   = "F1"
write  = false
nhits  = 128
Nf = 2

scalar  = false
deriv   = false
Nsmear  = ("0","40","80")
Ns_conn = [(2,1)]
xlim    = (0,12.5)
plot_sing, plot_conn = true, false
plt3 =  plot_from_logs(fileCONN,fileDISC,h5file,name,Nsmear,nhits;write,scalar,Nf,xlim,Ns_conn,plot_sing,plot_conn,deriv)
plot!(plt3;xlims=xlim,ylims=(0.50,0.70),legend=:topleft)

scalar  = true
deriv   = false
Nsmear  = ("0","40","80")
Ns_conn = []
xlim    = (0,6.5)
plot_sing, plot_conn = true, false
plt4 =  plot_from_logs(fileCONN,fileDISC,h5file,name,Nsmear,nhits;write,scalar,Nf,xlim,Ns_conn,plot_sing,plot_conn,deriv)
plot!(plt4;xlims=xlim,ylims=(0.50,0.80),legend=:topleft)

for plt in [plt3, plt4]
    add_mass_band!(plt,mπREF,ΔmπREF;label="PS",alpha=0.5)
    add_mass_band!(plt,mVREF,ΔmVREF;label="V",alpha=0.5)
    add_mass_band!(plt,mηREF,ΔmηREF;label="η",alpha=0.5)
    display(plt)
end

path = "/home/fabian/Documents/Physics/Data/DataDiaL/measurementsAS/Lt56Ls24beta6.8mas-1.035/out"
path = "/home/fabian/Dokumente/Physics/Data/DataDiaL/measurementsAS/Lt56Ls24beta6.8mas-1.035/out"
h5file = "test.hdf5"
mπREF, ΔmπREF = 0.3463, 0.0007
mρREF, ΔmρREF = 0.3989, 0.0016
mSREF, ΔmSREF = 0.5140, 0.0070
mηREF, ΔmηREF = NaN, NaN
write = false
nhits = 200
Nf = 3

fileCONN = joinpath(path,"out_spectrum_smeared")
fileDISC = joinpath(path,"out_spectrum_smeared_discon")
name     = "AS4"
Nsmear   = ("0","30","60")

fileCONN = joinpath(path,"out_spectrum_smeared_N240")
fileDISC = joinpath(path,"out_spectrum_smeared_discon_N240")
name     = "AS5"

scalar  = false
deriv   = !scalar 
Nsmear  = collect(0:240:240)
Ns_conn = [(2,1)]
xlim    = (0,12.5)
plot_sing, plot_conn = true, false
plt1 =  plot_from_logs(fileCONN,fileDISC,h5file,name,Nsmear,nhits;write,scalar,Nf,xlim,Ns_conn,plot_sing,plot_conn,deriv)
plot!(plt1;xlims=xlim,ylims=(0.30,0.45),legend=:topleft)

scalar  = true
deriv   = !scalar 
Nsmear  = (0,60,240)
Ns_conn = []
xlim    = (0,8.5)
plot_sing, plot_conn = true, false
plt2 =  plot_from_logs(fileCONN,fileDISC,h5file,name,Nsmear,nhits;write,scalar,Nf,xlim,Ns_conn,plot_sing,plot_conn,deriv)
plot!(plt2;xlims=xlim,ylims=(0.30,0.55),legend=:topleft)

for plt in [plt1,plt2]
    add_mass_band!(plt,mπREF,ΔmπREF;label="PS",alpha=0.5)
    add_mass_band!(plt,mρREF,ΔmρREF;label="V",alpha=0.5)
    add_mass_band!(plt,mSREF,ΔmSREF;label="S",alpha=0.5)
    display(plt)
end

#plot!(title="Effective mass: Scalar mesons (S) non-singlet")
#savefig("scalar_isovector.pdf")