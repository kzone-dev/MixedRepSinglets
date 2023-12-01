import gvar as gv
import corrfitter as cf
import h5py
import numpy as np

def get_hdf5_value(hdf5file,key):
    return hdf5file[key][()]

def make_models(T,tmin,tmax):
    """ Create corrfitter model for G(t). """
    return [cf.Corr2(datatag='Gab', tp=T, tmin=tmin, tmax=tmax, a='a', b='a', dE='dE')]

def make_prior(N):
    prior = gv.BufferDict()
    # NOTE: use a log-Gaussion distrubtion for forcing positive energies
    # NOTE: Even with this code they can be recovered by providing loose priors of 0.1(1) for both
    prior['log(a)']  = gv.log(gv.gvar(N * ['1(1)']))
    prior['log(dE)'] = gv.log(gv.gvar(N * ['1(1)']))
    return prior

def bootstrap_fit(fitter,dset,T,tmin,tmax,n=20,printing=False):
    pdatalist = (cf.process_dataset(ds, make_models(T,tmin,tmax)) for ds in gv.dataset.bootstrap_iter(dset, n=n))
    bs = gv.dataset.Dataset()
    for bsfit in fitter.bootstrapped_fit_iter(pdatalist=pdatalist):
        bs.append(E=np.cumsum(bsfit.pmean['dE']),a=bsfit.pmean['a'])
    bs = gv.dataset.avg_data(bs, bstrap=True)
    E = bs['E']
    a = bs['a']
    if printing:
        print('bootstrap: ',30 * '=')
        print('{:2}  {:15}  {:15}'.format('E', E[0], E[1]))
        print('{:2}  {:15}  {:15}'.format('a', a[0], a[1]))
    return E, a

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

def fit_correlator_with_bootstrap(data,T,tmin,tmax,Nmax,plotting=False,printing=False):
    T = abs(T) 
    fitter = cf.CorrFitter(models=make_models(T,tmin,tmax))
    avg = gv.dataset.avg_data(data)
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
    # NOTE: A bootstrap fit can only be performed if`the object `fitter` has 
    # already been used to perform a fit.
    # NOTE: The bootstrap analysis is performed using the priors and initial 
    # parameters used in the last invokation of the previous fit. 
    E_bs, a_bs = bootstrap_fit(fitter, data, T, tmin, tmax)
    # NOTE: From the lsqfit documentation
    # There are several different views available for each plot, specified by parameter view:os.
    #   'ratio': Data divided by fit (default).
    #   'diff': Data minus fit, divided by data’s standard deviation.
    #   'std': Data and fit.
    #   'log': 'std' with log scale on the vertical axis.
    #   'loglog': ‘std’` with log scale on both axes.
    if plotting:
        fit.show_plots(view='ratio')
        fit.show_plots(view='log'  )
    return E, a, E_bs, a_bs, chi2, dof

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

def fit_from_file(hdf5file,channel,tmin=0,delta_tmax=0,Nmax=10):
    f = h5py.File(hdf5file)
    T    = get_hdf5_value(f,"lattice")[0]
    L    = get_hdf5_value(f,"lattice")[1]
    corr = get_hdf5_value(f,channel)
    beta = get_hdf5_value(f,"beta")
    mass = get_hdf5_value(f,"quarkmasses")[0]

    dset = gv.dataset.Dataset(dict(Gab=np.transpose(corr)))
    # include binning
    dset = gv.dataset.bin_data(dset,binsize=2)
    # renormalization from lattice perturbation theory 
    tmax = T/2 - delta_tmax    
    E, a, E_bs, a_bs, chi2, dof = fit_correlator_with_bootstrap(dset,T,tmin,tmax,Nmax,plotting=True,printing=True)
    
    print("T,L,m0,beta,m_meson,chi2/dof")
    print( T,",",L,",",mass,",",beta,",",E[0],",",E[1],",",E[2],",",chi2/dof)

def fit_eigenvalues_file(hdf5file,tmin=0,delta_tmax=10,Nmax=10):
    f = h5py.File(hdf5file)
    T    = get_hdf5_value(f,"lattice")[0]
    L    = get_hdf5_value(f,"lattice")[1]
    ev = get_hdf5_value(f,"singlet_eigenvalues_g5")[()]
    Delta_ev = get_hdf5_value(f,"Delta_singlet_eigenvalues_g5")[()]
    beta = get_hdf5_value(f,"beta")
    mass = get_hdf5_value(f,"quarkmasses")[0]

    eig1 = dict(Gab=gv.gvar(ev[:,0],Delta_ev[:,0]))
    eig2 = dict(Gab=gv.gvar(ev[:,1],Delta_ev[:,1]))

    # renormalization from lattice perturbation theory 
    tmax = T/2 - delta_tmax    
    E1, a1, chi2A, dofA = fit_correlator_without_bootstrap(eig1,T,tmin,tmax,Nmax,plotting=True,printing=True)
    E2, a2, chi2B, dofB = fit_correlator_without_bootstrap(eig2,T,tmin,tmax,Nmax,plotting=True,printing=True)
    
    print("T,L,m0,beta,m_meson,chi2/dof")
    print( T,",",L,",",mass,",",beta,",",E1[0],",",chi2A/dofA)
    print( T,",",L,",",mass,",",beta,",",E2[0],",",chi2B/dofB)


#infile = "/home/fabian/Documents/Analysis/MixedRepSinglets/output/h5files/Lt64Ls20beta6.5mf0.71mas1.01FUN/out_spectrum.h5"
#fit_from_file(infile,"g5",tmin=2,delta_tmax=2,Nmax=10)

infile = "/home/fabian/Documents/Analysis/MixedRepSinglets/output/eigenvalues/eigenvalues_Lt48Ls20beta6.5mf0.71mas1.01.h5"
fit_eigenvalues_file(infile,tmin=0,delta_tmax=0,Nmax=10)
