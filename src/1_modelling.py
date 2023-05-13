'''
Stellar modelling with free values for surface correction.
'''
import os
import numpy as np
import pandas as pd
import sys

from grid import grid, get_model_Dnu
from astropy.io import ascii
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
# import numpy.lib.recfunctions
import h5py
import scipy
import scipy.interpolate

def read_models(filepath):
    atrack = Table(np.load(filepath, allow_pickle=True))
    return atrack

if __name__ == '__main__':

    ### step 0: gather the input parameters
    ipart = int(sys.argv[1]) # the no. of chunk processed
    Nthread = 8 # the number of parallel processing threads
    work_dir = './' # input data are stored here
    suffix = 'spec_gaia_nu' # output results are stored here
    model = 'combined' # the surface correction model, 'combined' or 'cubic'
    Nstar = 320 # number of stars processed at one time
    offset = 0
    outdir = '/Users/yaguang/results_{:s}/'.format(suffix) # the place to store results

    print('Processing stars from {} to {}'.format(ipart*Nstar+offset, (ipart+1)*Nstar+offset))
    if not os.path.exists(outdir): os.makedirs(outdir, exist_ok=True)
    random_index = int(np.random.rand()*1e8)

    
    ### step 1: prepare samples
    stars = pd.read_excel(work_dir+'sample/samples.xlsx')
    # kics = np.loadtxt('exists.txt', dtype='int')
    stars['luminosity'] = stars['lum_J']
    stars['e_luminosity'] = stars['e_lum_J']
    idx = ~np.isfinite(stars['Dnu'])
    stars.loc[idx, 'Dnu'] = stars.loc[idx, 'Dnu_guess']

    idx = np.isfinite(stars['ifmodelling']) & np.isfinite(stars['Dnu']) & (np.isfinite(stars['luminosity'])) & \
        (stars['[M/H]'] > -0.7) & np.isin(stars['names'], ['ngc6791', 'ngc6819'])  # & (~np.isin(stars['KIC'], kics)) 
    stars = stars.loc[idx,:].sort_values('[M/H]').reset_index(drop=True)
    idx = (stars.index>=(ipart*Nstar+offset)) & (stars.index<((ipart+1)*Nstar+offset))
    stars = stars.loc[idx,:].reset_index(drop=True)

    mh_llim, mh_ulim = np.min(stars['[M/H]'])-5*0.05, np.max(stars['[M/H]'])+5*0.05

    a, b, c = [1.3579, 1.2377, 3.2293]
    stars['mod_e_freq'] = a * (stars['numax']/3090)**b * (stars['Teff']/5777)**c

    cols = ['KIC', 'Teff', 'e_Teff', '[M/H]', 'e_[M/H]', 'luminosity', 'e_luminosity', 'Dnu', 'e_Dnu', 'numax', 'e_numax', 'mod_e_freq']
    stars.loc[:,cols].to_csv(work_dir+'inputs/samples_{:0.0f}.csv'.format(random_index), index=False)


    modes = pd.read_excel(work_dir+'sample/modes.xlsx')
    idx = (np.isin(modes['KIC'], stars['KIC'])) & (np.isin(modes['l'], [0])) # & 
    modes = modes.loc[idx,:].reset_index(drop=True)

    cols = ['KIC', 'l', 'fc', 'e_fc']
    modes.loc[:,cols].to_csv(work_dir+'inputs/modes_{:0.0f}.csv'.format(random_index), index=False)


    ### step 3: prepare a list of evolutionary tracks
    lists = []
    tracks = ['/Users/yaguang/uncalibrated/'+f for f in os.listdir('/Users/yaguang/uncalibrated') if (f.endswith('.npy') )]


    ### step 4: setup params
    from grid import user_setup_params

    new_params = {
    "observables": ['luminosity', 'Teff', '[M/H]' ] , #
    "estimators": ['index', 'star_age', 'star_mass', 'radius', 'density',\
                'Teff','luminosity','FeH', 'numax_scaling',\
                'fov_core', 'fov_shell', 'amlt','Yinit','Zinit',\
                'profile_number'],

    "filepath_stellar_params": work_dir+'inputs/samples_{:0.0f}.csv'.format(random_index), 
    "col_starIDs": 'KIC',
    "filepath_stellar_freqs": work_dir+'inputs/modes_{:0.0f}.csv'.format(random_index), 

    # estimator names to appear in corner plots. suggests only keep essentials.
    "estimators_to_plot": ['star_age', 'star_mass', 'radius', 'density'],
    # if True, output top models to an h5 file
    "if_data": True, 
    "Nthread": Nthread, 

    "if_plot": True, 
    "if_plot_echelle": True,

    # filepath to output results 
    "filepath_output": outdir,

    "if_classical": True,
    "weight_classical": 1, 
    "if_seismic": True,
    "weight_seismic": 3, 
    "if_reduce_seis_chi2": True, 
    "if_correct_surface": True, 
    # added 
    "require_negative_surface_correction": False,
    "require_absolute_surface_correction_increase_with_nu": False,
    # added
    "surface_correction_formula": model, 
    "if_add_model_error": True,
    "col_mode_n": "mode_n", 
    # 2 - set by the column in the filepath_stellar_params file
    # ``col_model_error'' must be set.
    "add_model_error_method": 2, 
    # used to evaluate the systematic uncertainties. the model rms frequency difference at below percentile 
    # "rescale_percentile": 10, 
    # # used to set the systematic uncertainties from the filepath_stellar_params file
    "col_model_error": "mod_e_freq",
    "col_mode_freq": "mode_freq_o",  
    }


    for key in new_params.keys():
        user_setup_params[key] = new_params[key]



    # initialize a grid modelling class
    g = grid(read_models, tracks, user_setup_params)

    g.run()
