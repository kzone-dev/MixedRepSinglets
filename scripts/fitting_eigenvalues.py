import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import matplotlib.pyplot as plt

def get_hdf5_value(hdf5file,key):
    return hdf5file[key][()]

def make_models(T,tmin,tmax):
    """ Create corrfitter model for G(t). """
    return [cf.Corr2(datatag='Gab', tp=-T, tmin=tmin, tmax=tmax, a='a', b='a', dE='dE')]

def make_prior(N):
    prior = gv.BufferDict()
    # NOTE: use a log-Gaussion distrubtion for forcing positive energies
    # NOTE: Even with this code they can be recovered by providing loose priors of 0.1(1) for both
    prior['log(a)']  = gv.log(gv.gvar(N * ['0.1(0.1)']))
    prior['log(dE)'] = gv.log(gv.gvar(N * ['0.5(0.3)']))
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

def fit_correlator_without_bootstrap(avg,T,tmin,tmax,Nmax,plotting=False,printing=False):
    T = abs(T) 
    fitter = cf.CorrFitter(models=make_models(T,tmin,tmax))
    p0 = None
    # TODO: find good Nmax
    for N in range(2,Nmax+1):
        prior = make_prior(N)
        fit = fitter.lsqfit(data=avg, prior=prior, p0=p0)
        p0 = fit.pmean

        if printing:
            print('nterm =', N, 30 * '=')
            #print(fit)
            print_fit_param(fit)

    E, a, chi2, dof = first_fit_parameters(fit) 
    if plotting:
        fit.show_plots(view='ratio')
        fit.show_plots(view='log'  )
    return E, a, chi2, dof

def fit_eigenvalues_file(hdf5file,tmin,tmax1,tmax2,Nmax=10):
    f = h5py.File(hdf5file)
    T = get_hdf5_value(f,"M4/lattice")[0]
    L = get_hdf5_value(f,"M4/lattice")[1]

    ev = get_hdf5_value(f,"M4/eigvals_g5_singlet")[()]
    Delta_ev = get_hdf5_value(f,"M4/Delta_eigvals_g5_singlet")[()]

    eig1 = dict(Gab=gv.gvar(ev[:,17],Delta_ev[:,17]))
    eig2 = dict(Gab=gv.gvar(ev[:,16],Delta_ev[:,16]))

    E1, a1, chi2A, dofA = fit_correlator_without_bootstrap(eig1,T,tmin,tmax1,Nmax,plotting=PLOT,printing=PRINT)
    E2, a2, chi2B, dofB = fit_correlator_without_bootstrap(eig2,T,tmin,tmax2,Nmax,plotting=PLOT,printing=PRINT)
    
    beta = get_hdf5_value(f,"M4/beta")
    mass = get_hdf5_value(f,"M4/quarkmasses_fundamental")[0]

    print("T,L,m0,beta,m_meson,chi2/dof")
    print( T,",",L,",",mass,",",beta,",",E1[0],",",chi2A/dofA)
    print( T,",",L,",",mass,",",beta,",",E2[0],",",chi2B/dofB)

PLOT=False
PRINT=False

filename="/home/fabian/Downloads/smeared_singlet_eigenvalues_M1234.hdf5"
fit_eigenvalues_file(filename,3,10,10,2)