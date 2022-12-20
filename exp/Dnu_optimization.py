import numpy as np 
import matplotlib.pyplot as plt
import h5py
import corner
import asteroseismology as se

data = list(np.load('data/stellar_models_for_Dnu_optimization.npy',allow_pickle=True))


Nstars = len(data)
# old model: use all models for each star
# def model(thetas):
#     a, b, c, d = thetas
#     chi2 = np.zeros(Nstars)
#     for istar in range(Nstars):
#         numax, mass, radius, Teff = [data[istar][col].reshape(-1,1) for col in ['numax', 'mass', 'radius', 'Teff']]
#         sig = a * mass**b * radius**c * (Teff/5777)**d
#         weights = 1./(sig*(2*np.pi)) * np.exp(-(data[istar]['mod_freq']-numax)**2.0/sig**2.0) # =(1/sig^2)
#         try:
#             Nmodels, Nmodes = data[istar]['mod_freq'].shape
#             X = np.concatenate([np.ones((Nmodes,1)), np.arange(Nmodes).reshape(-1,1)], axis=1)
#             w = np.array([np.diagflat(weights[imod,:]) for imod in range(Nmodels)])
#             y = data[istar]['mod_freq'][:,:].reshape(Nmodels, Nmodes, 1)

#             Dnu_mods = (np.linalg.inv(X.T @ w @ X) @ ((X.T @ w ) @ y))[:,1,0]

#             res = se.quantile(Dnu_mods.reshape(-1,1), [0.16, 0.50, 0.84], weights=data[istar]['prob'].reshape(-1,1)).reshape(-1)
#             Dnu_mod = res[1]
#             e_Dnu_mod = (res[2]-res[1])/2.

#             chi2[istar] = (Dnu_mods-data[istar]['Dnu_syd'])**2.0/(e_Dnu_mod**2.0 + data[istar]['e_Dnu_syd']**2.0)
#         except np.linalg.LinAlgError:
#             chi2[istar] = np.nan
#             break
#     return np.sum(chi2)

# new model: use the best model
best = np.array([np.argmax(data[istar]['prob']) for istar in range(Nstars)])
def model(thetas):
    a, b, c, d = thetas
    chi2 = np.zeros(Nstars)
    for istar in range(Nstars):
        
#         numax, mass, radius, Teff = [data[istar][col][best[istar]] for col in ['numax', 'mass', 'radius', 'Teff']]
#         sig = a * mass**b * radius**c * (Teff/5777)**d
#         weights = 1./(sig*(2*np.pi)) * np.exp(-(data[istar]['mod_freq'][best[istar]]-numax)**2.0/sig**2.0) # =(1/sig^2)
#         try:
#             _, Nmodes = data[istar]['mod_freq'].shape
#             X = np.concatenate([np.ones((Nmodes,1)), np.arange(Nmodes).reshape(-1,1)], axis=1)
#             w = np.diagflat(weights)
#             y = data[istar]['mod_freq'][best[istar],:]

#             Dnu_mods = (np.linalg.inv(X.T @ w @ X) @ ((X.T @ w @ y.T)))[1]

#             chi2[istar] = (Dnu_mods-data[istar]['Dnu_syd'])**2.0/(data[istar]['e_Dnu_syd']**2.0)
#         except np.linalg.LinAlgError:
#             chi2[istar] = np.nan
#             break
        try:
            numax, mass, radius, Teff = [data[istar][col] for col in ['numax', 'mass', 'radius', 'Teff']]
            sig = a * mass**b * radius**c * (Teff/5777)**d 
            sig = sig.reshape(-1,1)
            numax = numax.reshape(-1,1)
            weights = 1./(sig*(2*np.pi)**0.5) * np.exp(-(data[istar]['mod_freq']-numax)**2.0/(2*sig**2.0))  # (Nmodel, Nmodes)
            # weights = weights/np.sum(weights, axis=1).reshape(-1,1) # (Nmodel, Nmodes)
            Nmodels, Nmodes = data[istar]['mod_freq'].shape
            
            Ns = np.arange(Nmodes).reshape(1,-1)
            mod_freqs = data[istar]['mod_freq'][:,:] # (Nmodel, Nmodes)
            num = (-np.sum(weights*Ns, axis=1)*np.sum(weights*mod_freqs, axis=1) + \
             np.sum(weights, axis=1)*np.sum(weights*Ns*mod_freqs, axis=1) )
            den = (np.sum(weights, axis=1)*np.sum(weights*Ns**2., axis=1) - np.sum(weights*Ns, axis=1)**2. )

            Dnu_mod = np.average(num/den, weights=data[istar]['prob']) 
            chi2[istar] = (Dnu_mod-data[istar]['Dnu_syd'])**2.0/(data[istar]['e_Dnu_syd']**2.0)

        except np.linalg.LinAlgError:
            chi2[istar] = np.nan
            break
    return np.sum(chi2)

params_range = [[1,325], [-5, 5], [-5, 5], [-5, 5]]

def log_prob(thetas):
    for itheta, theta in enumerate(thetas):
        if not (params_range[itheta][0] < theta < params_range[itheta][1]):
            return -np.inf
    m = model(thetas)
    if np.isfinite(m): 
        return -model(thetas)/2.
    else:
        return -np.inf

rand = np.random.randn(20, 4)
rand[:,0] = rand[:,0]*20
rand[:,1] = rand[:,1]*2

coords = np.array([255,0.26,-1.67,1.47]) + rand
nwalkers, ndim = coords.shape

import emcee

if __name__ == '__main__':
    from multiprocessing import Pool
    with Pool(4) as pool:
        # Initialize the sampler
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)

        max_n = 1000
        
        sampler.run_mcmc(coords, max_n, progress=True)
        
        samples = sampler.get_chain(flat=True)
        np.save('data/MCMC_samples_for_Dnu_optimization', samples)