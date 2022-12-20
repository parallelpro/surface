import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
rootpath = os.getenv('WORK_DIR')+'numax-sc-metallicity/'
sys.path.append(rootpath)
work_dir = rootpath+'surface/'
import h5py
import asteroseismology as se
from astropy.io import ascii
from astropy.table import Table
import corner
import matplotlib.colors
import scipy
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

# from PTMCMCSampler import PTMCMCSampler