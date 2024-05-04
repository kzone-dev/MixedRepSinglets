import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import matplotlib.pyplot as plt
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
    # TODO: find good Nmax
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

def fit_eigenvalues(outfile,outfileHR,hdf5file,tmin1,tmax1,tp,Nmax,ensemble,channel,rep,header=False):
    f = h5py.File(hdf5file)
    T = get_hdf5_value(f,ensemble+"/"+rep+"/CONN/lattice")[0]
    L = get_hdf5_value(f,ensemble+"/"+rep+"/CONN/lattice")[1]
    tp = tp*T if tp != 0 else None

    h5name = ensemble+"/"+rep+"/CONN/DEFAULT_SEMWALL TRIPLET/"+channel
    print(h5name)
    corr = get_hdf5_value(f,h5name)[()]
    
    dset = gv.dataset.avg_data(np.transpose(corr))
    eig1 = dict(Gab=dset)
    
    E1, a1, chi2A, dofA = fit_correlator_without_bootstrap(eig1,T,tmin1,tmax1,Nmax,tp,plotting=PLOT,printing=PRINT)
    
    beta = get_hdf5_value(f,ensemble+"/"+rep+"/CONN/beta")

    out = open(outfile, "a")
    outHR = open(outfileHR, "a")
    out.write("%s;%s;%s;%s;%s;%s;%s;%s\n" % (ensemble,channel,T,L,beta,gv.mean(E1[0]),gv.sdev(E1[0]),chi2A/dofA))
    outHR.write("%s;%s;%s;%s;%s;%s;%s\n" % (ensemble,channel,T,L,beta,E1[0],chi2A/dofA))
    out.close()
    outHR.close()


def run_corrfitter_singlets(prmfile,hdf5file,outdir):
    outfile    = os.path.join(outdir,"corrfitter_results.csv")
    outfileHR  = os.path.join(outdir,"corrfitter_results_HR.csv")
    os.path.exists(outfile)   and os.remove(outfile)
    os.path.exists(outfileHR) and os.remove(outfileHR)
    

    with open(prmfile) as csvfile:
        reader = csv.DictReader(csvfile,delimiter=';')
        for row in reader:
            ensemble, channel, rep = row['ensemble'], row['channel'], row["rep"]
            tmin, tmax = int(row['tmin']), int(row['tmax'])
            tp,   Nmax = int(row['tp']), int(row['Nmax'])
            fit_eigenvalues(outfile,outfileHR,hdf5file,tmin,tmax,tp,Nmax,ensemble,channel,rep)


PLOT=True
PRINT=False

args = sys.argv
if len(args) < 4:
    print("Missing parameter and/or hdf5 file")
elif len(args)>=4:
    prmfile  = args[1]
    hdf5path = args[2] 
    outdir   = args[3] 
    run_corrfitter_singlets(prmfile,hdf5path,outdir)
