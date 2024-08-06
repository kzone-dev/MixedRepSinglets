import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import csv
import os
import sys

def get_hdf5_value(hdf5file,key):
    return hdf5file[key][()]

def make_models(T,tmin,tmax,tp):
    """ Create corrfitter model for G(t). """
    return [cf.Corr2(datatag='Gab', tp=tp, tmin=tmin, tmax=tmax, a='a', b='a', dE='dE')]

def make_prior(N):
    prior = gv.BufferDict()
    # setting the sdev of the prioir to infinity amounts to turning off the prior contribution to chi2
    prior['log(a)']  = gv.log(gv.gvar(N*[0.1], N*[0.1]))
    prior['log(dE)'] = gv.log(gv.gvar(N*[0.1], N*[0.1]))
    return prior

def first_fit_parameters(fit):
    p = fit.p
    E = np.cumsum(p['dE'])
    a = p['a']
    chi2 = fit.chi2     
    dof = fit.dof
    return E, a, chi2, dof

def fit_correlator_without_bootstrap(avg,T,tmin,tmax,Nmax,tp,plotting=False,printing=False):
    T = abs(T) 
    fitter = cf.CorrFitter(models=make_models(T,tmin,tmax,tp))
    p0 = None
    for N in range(1,Nmax+1):
        prior = make_prior(N)
        fit = fitter.lsqfit(data=avg, prior=prior, p0=p0)
        p0 = fit.pmean

        if printing:
            print('nterm =', N, 30 * '=')
            print(fit)

    E, a, chi2, dof = first_fit_parameters(fit) 
    if plotting:
        fit.show_plots(view='ratio')
        fit.show_plots(view='log'  )
    return E, a, chi2, dof

def fit_eigenvalues(outfile,outfileHR,hdf5file,tmin1,tmin2,tmax1,tmax2,tp,Nmax,ensemble,channel,header=False):
    f = h5py.File(hdf5file)
    T = get_hdf5_value(f,ensemble+"/"+channel+"/lattice")[0]
    L = get_hdf5_value(f,ensemble+"/"+channel+"/lattice")[1]
    tp = tp*T if tp != 0 else None

    ev = get_hdf5_value(f,ensemble+"/"+channel+"/eigvals")[()]
    cov_ev = get_hdf5_value(f,ensemble+"/"+channel+"/eigvals_cov")[()]
    Delta_ev = get_hdf5_value(f,ensemble+"/"+channel+"/Delta_eigvals")[()]

    Nops = ev.shape[1]
    eig1 = dict(Gab=gv.gvar(ev[:,Nops-1],cov_ev[:,:,Nops-1]))
    eig2 = dict(Gab=gv.gvar(ev[:,Nops-2],cov_ev[:,:,Nops-2]))

    E1, a1, chi2A, dofA = fit_correlator_without_bootstrap(eig1,T,tmin1,tmax1,Nmax,tp,plotting=PLOT,printing=PRINT)
    E2, a2, chi2B, dofB = fit_correlator_without_bootstrap(eig2,T,tmin2,tmax2,Nmax,tp,plotting=PLOT,printing=PRINT)
    
    beta = get_hdf5_value(f,ensemble+"/"+channel+"/beta")
    mf   = get_hdf5_value(f,ensemble+"/"+channel+"/quarkmasses_fundamental")[0]
    mas  = get_hdf5_value(f,ensemble+"/"+channel+"/quarkmasses_antisymmetric")[0]

    out = open(outfile, "a")
    outHR = open(outfileHR, "a")
    out.write("%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n" % (ensemble,channel,T,L,mf,mas,beta,gv.mean(E1[0]),gv.sdev(E1[0]),gv.mean(E2[0]),gv.sdev(E2[0]),chi2A/dofA,chi2B/dofB))
    outHR.write("%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n" % (ensemble,channel,T,L,mf,mas,beta,E1[0],E2[0],chi2A/dofA,chi2B/dofB))
    out.close()
    outHR.close()

def fit_eigenvalues_resample(outfile,outfileHR,hdf5file,tmin1,tmin2,tmax1,tmax2,tp,Nmax,ensemble,channel,header=False):
    f = h5py.File(hdf5file)
    T = get_hdf5_value(f,ensemble+"/"+channel+"/lattice")[0]
    L = get_hdf5_value(f,ensemble+"/"+channel+"/lattice")[1]
    tp = tp*T if tp != 0 else None

    Delta_ev = get_hdf5_value(f,ensemble+"/"+channel+"/Delta_eigvals")[()]
    cov_ev   = get_hdf5_value(f,ensemble+"/"+channel+"/eigvals_cov")[()]
    ev_resamples = get_hdf5_value(f,ensemble+"/"+channel+"/eigvals_resamples")[()]
    
    T, nsamples, Nops = ev_resamples.shape

    E1_samples = np.zeros(nsamples)
    E2_samples = np.zeros(nsamples)
    a1_samples = np.zeros(nsamples)
    a2_samples = np.zeros(nsamples)
    chi2dofA_samples = np.zeros(nsamples)
    chi2dofB_samples = np.zeros(nsamples)

    for n in range(nsamples):

        eig1 = dict(Gab=gv.gvar(ev_resamples[:,n,Nops-1],cov_ev[:,:,Nops-1]))
        eig2 = dict(Gab=gv.gvar(ev_resamples[:,n,Nops-2],cov_ev[:,:,Nops-2]))

        E1, a1, chi2A, dofA = fit_correlator_without_bootstrap(eig1,T,tmin1,tmax1,Nmax,tp,plotting=False,printing=False)
        E2, a2, chi2B, dofB = fit_correlator_without_bootstrap(eig2,T,tmin2,tmax2,Nmax,tp,plotting=False,printing=False)

        E1_samples[n] = gv.mean(E1[0]) 
        E2_samples[n] = gv.mean(E2[0])
        a1_samples[n] = gv.mean(a1[0])
        a2_samples[n] = gv.mean(a2[0])
        chi2dofA_samples[n] = chi2A/dofA 
        chi2dofB_samples[n] = chi2B/dofB

    # this dataset uses jackknife resampling
    E1 = gv.gvar(np.mean(E1_samples), np.sqrt(nsamples)*np.std(E1_samples,ddof=1))
    E2 = gv.gvar(np.mean(E2_samples), np.sqrt(nsamples)*np.std(E2_samples,ddof=1))
    a1 = gv.gvar(np.mean(a1_samples), np.sqrt(nsamples)*np.std(a1_samples,ddof=1))
    a2 = gv.gvar(np.mean(a2_samples), np.sqrt(nsamples)*np.std(a2_samples,ddof=1))
    chi2dofA = np.mean(chi2dofA_samples)
    chi2dofB = np.mean(chi2dofB_samples)

    beta = get_hdf5_value(f,ensemble+"/"+channel+"/beta")
    mf   = get_hdf5_value(f,ensemble+"/"+channel+"/quarkmasses_fundamental")[0]
    mas  = get_hdf5_value(f,ensemble+"/"+channel+"/quarkmasses_antisymmetric")[0]

    out = open(outfile, "a")
    outHR = open(outfileHR, "a")
    out.write("%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n" % (ensemble,channel,T,L,mf,mas,beta,gv.mean(E1),gv.sdev(E1),gv.mean(E2),gv.sdev(E2),chi2dofA,chi2dofB))
    outHR.write("%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n" % (ensemble,channel,T,L,mf,mas,beta,E1,E2,chi2dofA,chi2dofB))
    out.close()
    outHR.close()

def run_corrfitter_singlets(prmfile,hdf5path,outdir,resample=False):

    outfile    = os.path.join(outdir,"corrfitter_results.csv")
    outfileHR  = os.path.join(outdir,"corrfitter_results_HR.csv")
    outfile2   = os.path.join(outdir,"corrfitter_results_jackknife.csv")
    outfile2HR = os.path.join(outdir,"corrfitter_results_jackknife_HR.csv")
    os.path.exists(outfile)     and os.remove(outfile)
    os.path.exists(outfileHR)   and os.remove(outfileHR)
    os.path.exists(outfile2)    and os.remove(outfile2)
    os.path.exists(outfile2HR)  and os.remove(outfile2HR)
    
    hdf5file = os.path.join(hdf5path,"singlets_smeared_eigenvalues.hdf5")

    with open(prmfile) as csvfile:
        reader = csv.DictReader(csvfile,delimiter=';')
        for row in reader:
            ensemble, channel = row['ensemble'], row['channel']
            tmin1, tmin2 = int(row['tmin1']), int(row['tmin2'])
            tmax1, tmax2 = int(row['tmax1']), int(int(row['tmax2']))
            tp, Nmax = int(row['tp']), int(row['Nmax'])

            fit_eigenvalues(outfile,outfileHR,hdf5file,tmin1,tmin2,tmax1,tmax2,tp,Nmax,ensemble,channel)
            if resample:
                fit_eigenvalues_resample(outfile2,outfile2HR,hdf5file,tmin1,tmin2,tmax1,tmax2,tp,Nmax,ensemble,channel)

PLOT=False
PRINT=False

args = sys.argv
if len(args) < 4:
    print("Missing parameter and/or hdf5 file")
elif len(args)>=5:
    prmfile  = args[1]
    hdf5path = args[2] 
    outdir   = args[3] 
    resample = args[4] == 'True'
    run_corrfitter_singlets(prmfile,hdf5path,outdir,resample)
elif len(args)==4:
    prmfile  = args[1]
    hdf5path = args[2] 
    outdir   = args[3] 
    run_corrfitter_singlets(prmfile,hdf5path,outdir,resample=False)
