using HiRepParsing
using MixedRepSinglets
using LaTeXStrings
using Plots
include("utils.jl")     
gr(frame=:box,legendfontsize=8,legend=:bottomright)
plotlyjs(frame=:box,legendfontsize=8,legend=:bottomright)
pgfplotsx(legendfonthalign = :left,legend=:topright, frame=:box, legendfontsize=12, tickfontsize=14, labelfontsize=14, titlefontsize=16,  markersize=5)

function correlation_matrices_single_rep(fileCONN,fileDISC,h5file,name,Nsmear,nhits;write,scalar,Nf,deriv)

    if write
        typesCONN = r"source_N[0-9]+_sink_N[0-9]+ TRIPLET"
        typesDISC = r"DISCON_SEMWALL smear_N[0-9]+ SINGLET"
        writehdf5_spectrum_with_regexp(fileCONN,h5file,typesCONN,h5group="$name/CONN")
        writehdf5_spectrum_disconnected_with_regexp(fileDISC,h5file,typesDISC,nhits,h5group="$name/DISC")
    end

    if scalar
        conn, disc = _assemble_correlation_matrix_rep(h5file,name,Nsmear,"";channel="id",disc_sign=-1,subtract_vev=true, Nf,nsrc_max=nhits÷1)
    else
        conn, disc = _assemble_correlation_matrix_rep(h5file,name,Nsmear,"";channel="g5",disc_sign=+1,subtract_vev=false,Nf,nsrc_max=nhits÷1)
    end
    conn = correlator_folding(conn;t_dim=4,sign=+1)
    disc = correlator_folding(disc;t_dim=4,sign=+1)
    correlation_matrix = @. conn - disc
    correlation_matrix_deriv = correlator_derivative(correlation_matrix,t_dim=4)
    
    binsize = 1
    t0      = 2 # at least 2 needed for scalar states
    eigvals1, Δeigvals1, meff1, Δmeff1, eigenvalues_jackknife1 = eigenvalues_meff_mixed_rep(correlation_matrix;t0,binsize,deriv)
    eigvals2, Δeigvals2, meff2, Δmeff2, eigenvalues_jackknife2 = eigenvalues_meff_mixed_rep(conn;t0,binsize,deriv=false)
    return meff1, Δmeff1, meff2, Δmeff2, conn, correlation_matrix, correlation_matrix_deriv 
end
function plot_singlet_vs_nonsinglets(meff1, Δmeff1, meff2, Δmeff2, conn, correlation_matrix, correlation_matrix_deriv, Nsmear; xlim, Ns_conn, deriv, plot_sing=false, plot_conn=false)

    Nl  = length(Nsmear)
    plt = plot(ylabel=L"m_{eff}",xlabel=L"t")
    shapes = [:pentagon, :rect, :hexagon]

    t = 1:Int(round(xlim[2]))

    plot_sing && plot!(plt,t,meff1[Nl,t], yerr = Δmeff1[Nl,t],label="singlet (GEVP)",ms=5,markershape=:circle)
    plot_conn && plot!(plt,t,meff2[Nl,t], yerr = Δmeff2[Nl,t],label="non-singlet (GEVP)",ms=5,markershape=:rect)

    for (i,ind) in enumerate(Ns_conn) 
        corr  = deriv ? correlation_matrix_deriv[ind[1],ind[2],:,:] : correlation_matrix[ind[1],ind[2],:,:]
        corr0 = conn[ind[1],ind[2],:,:]
        sign  = deriv ? -1 : +1
        
        meff , Δmeff  = implicit_meff_jackknife(corr' ;sign)
        meff0, Δmeff0 = implicit_meff_jackknife(corr0';sign=+1)
        smear="(N1=$(Nsmear[ind[1]]), N2=$(Nsmear[ind[2]]))"

        plot_sing && plot!(plt,t,meff[t] , yerr = Δmeff[t],  ms=5, markershape=shapes[i]; label="singlet: $smear")
        plot_conn && plot!(plt,t,meff0[t], yerr = Δmeff0[t], ms=5, markershape=shapes[i]; label="non-singlet: $smear")
    end

    plot!(plt;xlims=xlim,ylims=(0.30,0.45),legend=:topleft)
    return plt
end

path = "/home/fabian/Downloads/smearing_old_complex_numbers/"
fileCONN = joinpath(path,"out_spectrum_smeared")
fileDISC = joinpath(path,"out_spectrum_smeared_discon")
name  = "F1"
nhits = 128
h5file = "test.hdf5"
Nsmear = ("0","40")
Nsmear = ("0","40","80")
mπREF, ΔmπREF = 0.5663, 0.0008  
mηREF, ΔmηREF = 0.6100, 0.0060  
Nf = 2

path = "/home/fabian/Documents/Physics/Data/DataDiaL/measurementsAS/Lt56Ls24beta6.8mas-1.035/out"
path = "/home/fabian/Dokumente/Physics/Data/DataDiaL/measurementsAS/Lt56Ls24beta6.8mas-1.035/out"
h5file = "test.hdf5"
mπREF, ΔmπREF = 0.3463, 0.0007
mρREF, ΔmρREF = 0.3989, 0.0016
mSREF, ΔmSREF = 0.5140, 0.0070
mηREF, ΔmηREF = NaN, NaN
nhits = 200
Nf = 3

fileCONN = joinpath(path,"out_spectrum_smeared")
fileDISC = joinpath(path,"out_spectrum_smeared_discon")
name     = "AS4"
Nsmear   = ("0","30","60")

fileCONN = joinpath(path,"out_spectrum_smeared_N240")
fileDISC = joinpath(path,"out_spectrum_smeared_discon_N240")
name     = "AS5"

Nsmear   = collect(0:240:240)
Ns_conn  = [(2,1)]
xlim     = (0,12.5)

write  = false
scalar = false
deriv  = true 
plot_sing, plot_conn = true, false

meff1, Δmeff1, meff2, Δmeff2, conn, correlation_matrix, correlation_matrix_deriv = correlation_matrices_single_rep(fileCONN,fileDISC,h5file,name,Nsmear,nhits;write,scalar,Nf,deriv)
plt = plot_singlet_vs_nonsinglets(meff1, Δmeff1, meff2, Δmeff2, conn, correlation_matrix, correlation_matrix_deriv, Nsmear; xlim, Ns_conn, deriv, plot_sing, plot_conn)

add_mass_band!(plt,mπREF,ΔmπREF;label="PS",alpha=0.5)
add_mass_band!(plt,mρREF,ΔmρREF;label="V",alpha=0.5)
#add_mass_band!(plt,mSREF,ΔmSREF;label="S",alpha=0.5)
#add_mass_band!(plt,mηREF,ΔmηREF;label="",alpha=0.5)

display(plt)
#plot!(title="Effective mass: Scalar mesons (S) non-singlet")
#savefig("scalar_isovector.pdf")