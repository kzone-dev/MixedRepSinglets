import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import matplotlib.pyplot as plt
import csv

def get_hdf5_value(hdf5file,key):
    return hdf5file[key][()]

def make_models(T,tmin,tmax,tp):
    """ Create corrfitter model for G(t). """
    return [cf.Corr2(datatag='Gab', tp=tp, tmin=tmin, tmax=tmax, a='a', b='a', dE='dE')]

def make_prior(N):
    prior = gv.BufferDict()
    # NOTE: use a log-Gaussion distrubtion for forcing positive energies
    # NOTE: Even with this code they can be recovered by providing loose priors of 0.1(1) for both
    prior['log(a)']  = gv.log(gv.gvar(N * ['1(1)']))
    prior['log(dE)'] = gv.log(gv.gvar(N * ['1(1)']))
    return prior

def first_fit_parameters(fit):
    p = fit.p
    E = np.cumsum(p['dE'])
    a = p['a']
    chi2 = fit.chi2     
    dof = fit.dof
    return E, a, chi2, dof

def print_fit_param(fit):
    E, a, chi2, dof = first_fit_parameters(fit) 
    print('{:2}  {:15}  {:15}'.format('E', E[0], E[1]))
    print('{:2}  {:15}  {:15}'.format('a', a[0], a[1]))
    print('chi2/dof = ', chi2/dof, '\n')

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
            print_fit_param(fit)

    E, a, chi2, dof = first_fit_parameters(fit) 
    if plotting:
        fit.show_plots(view='ratio')
        fit.show_plots(view='log'  )
    return E, a, chi2, dof

def fit_eigenvalues_file(hdf5file,tmin1,tmin2,tmax1,tmax2,tp,Nmax,ensemble,channel,header=False):
    f = h5py.File(hdf5file)
    T = get_hdf5_value(f,ensemble+"/"+channel+"/lattice")[0]
    L = get_hdf5_value(f,ensemble+"/"+channel+"/lattice")[1]
    tp = tp*T if tp != 0 else None

    ev = get_hdf5_value(f,ensemble+"/"+channel+"/eigvals")[()]
    Delta_ev = get_hdf5_value(f,ensemble+"/"+channel+"/Delta_eigvals")[()]

    Nops = ev.shape[1]
    eig1 = dict(Gab=gv.gvar(ev[:,Nops-1],Delta_ev[:,Nops-1]))
    eig2 = dict(Gab=gv.gvar(ev[:,Nops-2],Delta_ev[:,Nops-2]))

    E1, a1, chi2A, dofA = fit_correlator_without_bootstrap(eig1,T,tmin1,tmax1,Nmax,tp,plotting=PLOT,printing=PRINT)
    E2, a2, chi2B, dofB = fit_correlator_without_bootstrap(eig2,T,tmin2,tmax2,Nmax,tp,plotting=PLOT,printing=PRINT)
    
    beta = get_hdf5_value(f,ensemble+"/"+channel+"/beta")
    mass = get_hdf5_value(f,ensemble+"/"+channel+"/quarkmasses_fundamental")[0]

    if header: 
        print("ensemble,channel,T,L,m0,beta,m_meson,chi2/dof")
    
    print(ensemble,",",channel,",",T,",",L,",",mass,",",beta,",",E1[0],",",chi2A/dofA)
    print(ensemble,",",channel,",",T,",",L,",",mass,",",beta,",",E2[0],",",chi2B/dofB)

PLOT=False
PRINT=False

filename="/home/fabian/Downloads/smeared_singlet_eigenvalues_M1234_with_conn_v2.hdf5"

with open('input/parameters_corrfitter.csv') as csvfile:
    reader = csv.DictReader(csvfile,delimiter=';')
    for row in reader:
        ensemble, channel = row['ensemble'], row['channel']
        tmin1, tmin2 = int(row['tmin1']), int(row['tmin2'])
        tmax1, tmax2 = int(row['tmax1']), int(int(row['tmax2']))
        tp, Nmax = int(row['tp']), int(row['Nmax'])
        fit_eigenvalues_file(filename,tmin1,tmin2,tmax1,tmax2,tp,Nmax,ensemble,channel,header=False)
