using Pkg; Pkg.activate(".")
using MixedRepSinglets
using Plots
using HDF5
#gr(fontfamily="Computer Modern",  top_margin=4Plots.mm, left_margin=4Plots.mm, legend=:topright, frame=:box, legendfontsize=11, tickfontsize=10, labelfontsize=14, markersize=6)

h5file = "/home/fabian/Downloads/single_rep_smeared.hdf5"
h5data = "/home/fabian/Downloads/data.hdf5"

N_max  = 80
Nsmear = 0:40:N_max
T,L    = h5read(h5file,"FUN/CONN/lattice")[1:2]
Ns     = length(Nsmear)
sign   = +1
Nf     = 2

# read data with and without APE smearing
conn_matrix, disc_matrix = _assemble_correlation_matrix_rep(h5file,"",Nsmear,"FUN";channel="g5",Nf=Nf)

# add data from publication as reference
data_conn = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90/out_spectrum/DEFAULT_SEMWALL TRIPLET_g5"
data_disc = "runsSp4/Lt24Ls12beta6.9m1-0.90m2-0.90/out_spectrum_discon/DISCON_SEMWALL SINGLET_g5_disc_re"
conn_old = h5read(h5data,data_conn)
loop_old = h5read(h5data,data_disc)
disc_old = Nf*unbiased_estimator(loop_old;subtract_vev = false, rescale=L^3)
MixedRepSinglets.rescale_connected!(conn_old,L)

# built full singlet correlator
corr_new = conn_matrix - disc_matrix
corr_old = conn_old'   - disc_old

corr_new = correlator_folding(corr_new;t_dim=4,sign)
corr_old = correlator_folding(corr_old;t_dim=2,sign)

c1, Δc1 = stdmean(corr_new,dims=3)
c3, Δc3 = stdmean(corr_old,dims=1)

# symmetry sign of correlators
corr_new = correlator_derivative(corr_new;t_dim=4)
corr_old = correlator_derivative(corr_old;t_dim=2)
sign = -1

corr_new
corr_old

samples = eigenvalues_jackknife_samples(corr_new,t0=1)
meff_new, Δmeff_new = meff_from_jackknife(samples;sign)
meff_old, Δmeff_old = implicit_meff_jackknife(corr_old';sign)
meff_APE, Δmeff_APE = implicit_meff_jackknife(corr_new[1,1,:,:]';sign)

# Plot correlator data without drivative
plt1 = plot(yscale=:log10)
scatter!(plt1,c3       ,yerr=Δc3       , label ="N=0 (with APE)")
scatter!(plt1,c1[1,1,:],yerr=Δc1[1,1,:], label ="published data")
display(plt1)

# Plot effective mass
range  = 2:10 
xlabel = "Euclidean time"
ylabel = "effective mass"
title  = "pseudoscalar singlet: with numerical derivative"
plt2 = plot(;xlabel,ylabel,title)
plot!(plt2,xlims=(0,12),ylims=(0.5,0.9))
scatter!(plt2,range, meff_new[Ns,range],yerr=Δmeff_new[Ns,range], label="GEVP N=0...$N_max (with APE)")
scatter!(plt2,range, meff_old[range,],yerr=Δmeff_old[range,], label="no smearing")
scatter!(plt2,range, meff_APE[range,],yerr=Δmeff_APE[range,], label="only APE")
hspan!(plt2,[0.604,0.616],alpha=0.5,label="published result")
display(plt2)
