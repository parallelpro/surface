#! /usr/bin/env python

import numpy as np 
import pandas as pd 

class seismo():
    def __init__(self, path_grid_models='grid.parquet', 
                 Teff_sun=5777, numax_sun=3090, Dnu_sun=135.1, 
                 Z_X_sun=0.0181547, g_sun=27400, fnumax=1.,
                 dY_dZ=None, Yp=None, mean_Yinit=None, sig_Yinit=0.01,
                 mean_amlt=None, sig_amlt=0.1):
        '''
        Calculate stellar parameters (based on a grid of stellar models) given observables
        (among Teff, [M/H], Dnu, numax, lum and logg). The input physics of the stellar 
        models are described in Li, Y. et al. (2022, <ref>)
        
        Stellar model data reported in the paper can be downloaded from <URL>.
        The free input parameters are (mass, [M/H], Yinit, alpha_mixing_length).
        See the paper for a full description of the input physics.
        The parameter range is 0.8<M/Msun<2.2, -0.8<[M/H]<0.6, pre-RGB-tip Dnu/microHz>2.0
        
        Models applicable to metal-poor stars and later-evoltionary stages are under development.
        
        Usage:
        >>> observables = {'Teff':[5777,4671], 'e_Teff':[30,80], 
        >>>               'MH':[0.,0.29], 'e_MH':[0.05,0.15], 
        >>>               'numax':[3090,74.21], 'e_numax':[100,0.68], 
        >>>               'Dnu':[135.1,6.937], 'e_Dnu':[0.1,0.014],
        >>>               'lum':[1.0, np.nan], 'e_lum':[0.02, np.nan] }
        >>> s = seismo(path_grid_models='database/grid.parquet').set_observables(**observables)
        >>> print(s.get_mass())
        >>> print(s.get_scaling_mass())
        >>> print(s.get_radius())
        >>> print(s.get_scaling_radius())
        >>> print(s.get_age())
        
        '''
        
        # read in grid models
        self.grid = pd.read_parquet(path_grid_models)
    
        # set up solar reference values
        self.Teff_sun = Teff_sun
        self.numax_sun = numax_sun
        self.Dnu_sun = Dnu_sun 
        self.Z_X_sun = Z_X_sun 
        self.g_sun = g_sun
        
        self.Dnu_sun_model = 135.0959
        self.density_sun_model = 1.00570

        # set up fnumax correction factor 
        self.fnumax = fnumax
        
        # create columns according to the scaling relations 
        self.grid['lum'] = self.grid['radius']**2.0 * (self.grid['Teff']/self.Teff_sun)**4.0 
        self.grid['grav'] = self.grid['mass'] * self.grid['radius']**-2.
        self.grid['logg'] = np.log10(self.grid['grav']*self.g_sun)
        self.grid['density'] = self.grid['mass'] * self.grid['radius']**-3.
        self.grid['[M/H]'] = np.log10(self.grid['surface_Z/X']) - np.log10(self.Z_X_sun)
        self.grid['numax'] = self.grid['grav'] * (self.grid['Teff']/self.Teff_sun)**-0.5 * self.numax_sun
        
        self.grid['fDnu'] = (self.grid['Dnu_freq']/self.Dnu_sun_model) * (self.grid['density']/self.density_sun_model)**-0.5
        self.grid['fDnu_o'] = (self.grid['Dnu_freq_o']/self.Dnu_sun_model) * (self.grid['density']/self.density_sun_model)**-0.5
        # self.grid['fDnu_t'] = (self.grid['Dnu_freq_t']/self.Dnu_sun_model) * (self.grid['density']/self.density_sun_model)**-0.5
        
        # set up priors
        self.priors = np.ones(len(self.grid))
        
        # set priors on Yinit if needed, based on the chemical enrichment law or user-defined value
        if (dY_dZ is not None) and (Yp is not None) and (sig_Yinit is not None):
            Yinit_according_to_Z = Yp + dY_dZ * self.grid['Zinit'].to_numpy()
            Yinit = self.grid['Yinit'].to_numpy()
            self.priors *= np.exp(-((Yinit-Yinit_according_to_Z)**2.0/sig_Yinit**2.0)/2.)
        elif (mean_Yinit is not None) and (sig_Yinit is not None):
            Yinit = self.grid['Yinit'].to_numpy()
            self.priors *= np.exp(-((Yinit-mean_Yinit)**2.0/sig_Yinit**2.0)/2.)
        else:
            pass
        
        # set priors on amlt if needed, based on the user-defined value
        if (mean_amlt is not None) and (sig_amlt is not None):
            amlt = self.grid['amlt'].to_numpy()
            self.priors *= np.exp(-((amlt-mean_amlt)**2.0/sig_amlt**2.0)/2.)
        
        

        
    def set_observables(self, Teff=None, e_Teff=None, MH=None, e_MH=None,
                        Dnu=None, e_Dnu=None, numax=None, e_numax=None,
                        lum=None, e_lum=None, logg=None, e_logg=None,
                        mass=None, e_mass=None, radius=None, e_radius=None,
                        density=None, e_density=None):
        '''
        Pass the observables as 1d arrays. If the constraint for some stars is missing, then set it as np.nan
        '''
        cols = np.array(['Teff', '[M/H]', 'Dnu_freq', 'numax', 'lum', 'logg', 'mass', 'radius', 'density'])
        idx = np.array([not (var is None) for var in [Teff, MH, Dnu, numax, lum, logg, mass, radius, density]]) & \
              np.array([not (var is None) for var in [e_Teff, e_MH, e_Dnu, e_numax, e_lum, e_logg, e_mass, e_radius, e_density]])
        cols = cols[idx]
        self.Nobs = len(cols)
        
        obs_data = [np.atleast_1d(var) if (not (var is None)) else None for ivar, var in enumerate([Teff, MH, Dnu, numax, lum, logg, mass, radius, density]) ]
        e_obs_data = [np.atleast_1d(var) if (not (var is None)) else None for ivar, var in enumerate([e_Teff, e_MH, e_Dnu, e_numax, e_lum, e_logg, e_mass, e_radius, e_density]) ]
        self.Teff, self.MH, self.Dnu, self.numax, self.lum, self.logg, self.mass, self.radius, self.density = obs_data
        self.e_Teff, self.e_MH, self.e_Dnu, self.e_numax, self.e_lum, self.e_logg, self.e_mass, self.e_radius, self.e_density = e_obs_data

        obs_data = np.array([var for ivar, var in enumerate(obs_data) if idx[ivar] ]).reshape(1,self.Nobs,-1) # shape (1, Nobs, Nstar)
        e_obs_data = np.array([var for ivar, var in enumerate(e_obs_data) if idx[ivar] ]).reshape(1,self.Nobs,-1) # shape (1, Nobs, Nstar)
        
        mod_data = self.grid[cols].to_numpy().reshape(-1,self.Nobs,1) # shape (Nmod, Nobs, 1)
        self.probs = self.priors.reshape(-1,1) * np.exp(-(np.nansum((obs_data-mod_data)**2.0/e_obs_data**2.0, axis=1))/2.) # shape (Nmod, Nstar)
        
        idx = np.sum(self.probs, axis=1)>0
        self.probs = self.probs[idx, :]
        self.priors = self.priors[idx]
        self.grid = self.grid[idx]
        
        self.Nstar = obs_data.shape[2]
        
        return self
    
    
    def _quantile(self, x, q, weights):
        idx = np.argsort(x)
        sw = weights[idx]
        cdf = np.cumsum(sw)[:-1]
        cdf /= cdf[-1]
        cdf = np.append(0, cdf)
        return np.interp(q, cdf, x[idx])
    
    
    def _get_estimators(self, variable_name):
        samples = self.grid[variable_name].to_numpy()
        res = []
        for istar in range(self.Nstar):
            res.append(self._quantile(samples, (0.16, 0.5, 0.84), self.probs[:,istar]))
        res = np.array(res)
        
        return res[:,1], (res[:,2]-res[:,0])/2.
    
    def _out(self, *arrays):
        arrays = list(arrays)
        if self.Nstar==1:
            return [arr[0] for arr in arrays]
        else:
            return arrays
        
        
    def get_mass(self):
        '''
        Determine mass from the stellar models.
        '''
        return self._out(*self._get_estimators('mass'))
    
    
    def get_scaling_mass(self):
        '''
        Determine masses [Msun] with the seismic scaling relation (with fDnu calculated from the stellar models).
        Observables numax, Dnu, and Teff must be set.
        '''
        fDnu, e_fDnu = self._get_estimators('fDnu')
        scaling_mass = (self.numax/self.numax_sun)**3. * (self.fnumax)**-3. * (self.Dnu/self.Dnu_sun)**-4. * (fDnu)**4.0 * (self.Teff/self.Teff_sun)**1.5
        e_scaling_mass = scaling_mass * ((3*self.e_numax/self.numax)**2. + (4*self.e_Dnu/self.Dnu)**2. + (4*e_fDnu/fDnu)**2.0 + (1.5*self.e_Teff/self.Teff)**2.)**0.5
        return self._out(scaling_mass, e_scaling_mass)
    
    
    def get_radius(self):
        '''
        Determine radius [Rsun] from the stellar models.
        '''
        return self._out(*self._get_estimators('radius'))
    
    
    def get_scaling_radius(self):
        '''
        Determine radii [Rsun] with the seismic scaling relation (with fDnu calculated from the stellar models).
        Observables numax, Dnu, and Teff must be set.
        '''
        fDnu, e_fDnu = self._get_estimators('fDnu')
        scaling_radius = (self.numax/self.numax_sun)**1. * (self.fnumax)**-1. * (self.Dnu/self.Dnu_sun)**-2. * (fDnu)**2.0 * (self.Teff/self.Teff_sun)**0.5
        e_scaling_radius = scaling_radius * ((1*self.e_numax/self.numax)**2. + (2*self.e_Dnu/self.Dnu)**2. + (2*e_fDnu/fDnu)**2.0 + (0.5*self.e_Teff/self.Teff)**2.)**0.5
        return self._out(scaling_radius, e_scaling_radius)
    
    
    def get_fDnu(self):
        '''
        Determine fDnu from the stellar models.
        '''
        return self._out(*self._get_estimators('fDnu'))
    
    
    def get_grav(self):
        '''
        Determine grav [gsun] from the stellar models.
        '''
        return self._out(*self._get_estimators('grav'))
    
    
    def get_age(self):
        '''
        Determine age [Gyr] from the stellar models.
        '''
        return self._out(*self._get_estimators('age'))
    
    
    def get_logg(self):
        '''
        Determine logg [cgs] from the stellar models.
        '''
        return self._out(*self._get_estimators('logg'))
    
    
    def get_density(self):
        '''
        Determine grav [rhosun] from the stellar models.
        '''
        return self._out(*self._get_estimators('density'))

    
    def get_lum(self):
        '''
        Determine lum [Lsun] from the stellar models.
        '''
        return self._out(*self._get_estimators('lum'))
    
    
    def get_Dnu(self):
        '''
        Determine Dnu [muHz] from the stellar models (oscillation frequencies based, surface effect corrected).
        '''
        return self._out(*self._get_estimators('Dnu_freq'))
    
    
    def get_numax(self):
        '''
        Determine numax [muHz] from the stellar models (scaling relations based).
        '''
        return self._out(*self._get_estimators('numax'))
    
    

if __name__ == '__main__':
    
    observables = {'Teff':[5777,4671], 'e_Teff':[30,80], 
                   'MH':[0.,0.29], 'e_MH':[0.05,0.15], 
                   'numax':[3090,74.21], 'e_numax':[100,0.68], 
                   'Dnu':[135.1,6.937], 'e_Dnu':[0.1,0.014],
                   'lum':[1.0, np.nan], 'e_lum':[0.02, np.nan] }

    print('Case 0: the sun and an RGB star')
    s = seismo(path_grid_models='database/grid.parquet').set_observables(**observables)
    print(s.get_mass())
    print(s.get_scaling_mass())
    print(s.get_radius())
    print(s.get_scaling_radius())
    print(s.get_age())
    print(s.get_fDnu())
    print(s.get_grav())
    print(s.get_Dnu())
    print(s.get_density())
    print(s.get_numax())

    print('Case 1: another solar abundance scale')
    s = seismo(path_grid_models='database/grid.parquet', Z_X_sun=0.0225)
    s.set_observables(**observables)
    print(s.get_mass())
    print(s.get_scaling_mass())
    print(s.get_radius())
    print(s.get_scaling_radius())
    print(s.get_age())
    
    print('Case 2: applying chemical enrichment law as a prior')
    s = seismo(path_grid_models='database/grid.parquet', Yp=0.22, dY_dZ=2.1268681549909063, sig_Yinit=0.1)
    s.set_observables(**observables)
    print(s.get_mass())
    print(s.get_scaling_mass())
    print(s.get_radius())
    print(s.get_scaling_radius())
    print(s.get_age())
    