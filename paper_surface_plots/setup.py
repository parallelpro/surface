from matplotlib import rcParams
rcParams["figure.dpi"] = 100
rcParams["savefig.dpi"] = 100

import numpy as np
import corner
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from lightkurve import search_lightcurvefile
import lightkurve as lk
# import seaborn as sns
import os 
from matplotlib.colors import ListedColormap
import sys
from astropy.io import ascii
import asteroseismology as se 
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

import h5py
os.environ['TEXINPUTS'] = ".:$HOME/texmf/custom-styles//:"

# color stuff
# red = sns.xkcd_rgb["pale red"]
# blue = sns.xkcd_rgb["denim blue"]
# green = sns.xkcd_rgb["faded green"]
# orange = sns.xkcd_rgb["amber"]
# grey = sns.xkcd_rgb["greyish"]
# darkgrey = sns.xkcd_rgb["dark grey"]
# purple = sns.xkcd_rgb["dusty purple"]
# black = sns.xkcd_rgb["black"]

red = 'indianred' #sns.xkcd_rgb["pale red"]
blue = 'steelblue' #sns.xkcd_rgb["denim blue"]
lightblue = 'lightskyblue' #sns.xkcd_rgb["light blue"]
green = 'limegreen' #sns.xkcd_rgb["faded green"]
orange = 'darkorange' #sns.xkcd_rgb["amber"]
cyan = 'cyan' #sns.xkcd_rgb["cyan"]
grey = 'lightgray' #sns.xkcd_rgb["greyish"]
darkgrey = 'darkgray' #sns.xkcd_rgb["dark grey"]
purple = 'darkmagenta' #sns.xkcd_rgb["dusty purple"]
black = 'k' #sns.xkcd_rgb["black"]

# def cmap_diverging(n=10):
#     return ListedColormap(sns.diverging_palette(220, 20, n=n))

# def cmap_grey(n=10):
#     return ListedColormap(sns.color_palette("Greys", n))

# def blues(n=10):
#     return sns.color_palette("Blues", n)


rootpath = os.getenv('WORK_DIR')+'numax-sc-metallicity/'

def to_overleaf(figure_file_name, subdir=''):
     return "rclone copy {:s} remote:Apps/Overleaf/Yaguang_surface/{:s}/".format(figure_file_name, subdir)

    
overleaf_path = '/Users/yaguang/Dropbox (Personal)/Apps/Overleaf/Yaguang_surface/figures/'
work_path = rootpath+'Onedrive/Work/numax-sc-metallicity/'
sys.path.append(work_path)
# from lib.toolkit import *

fontsize = 7 # minimum fontsize is 5
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["font.size"] = fontsize #7.5
matplotlib.rcParams["legend.fontsize"] = fontsize#7.5
matplotlib.rcParams['text.usetex'] = True#False #True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{helvet}\renewcommand{\familydefault}{\sfdefault}\usepackage{sfmath}'
#plt.rc('font', family='serif')
matplotlib.rcParams['axes.labelsize'] = fontsize#7
matplotlib.rcParams['xtick.labelsize'] = fontsize#7
matplotlib.rcParams['ytick.labelsize'] = fontsize#7
matplotlib.rcParams['ytick.direction']='out'
matplotlib.rcParams['ytick.major.size']=3.0
matplotlib.rcParams['ytick.minor.size']=2.0
matplotlib.rcParams['xtick.direction']='out'
matplotlib.rcParams['xtick.major.size']=3.0
matplotlib.rcParams['xtick.minor.size']=2.0
matplotlib.rcParams['font.family'] = 'sans-serif' #'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica'] #'Helvetica' 




# mnras size in pt
columnwidth = 240
textwidth = 504

def mnras_size(column="one", square=False, ratio=None):
    # Thanks Dan!
    # Parameters:
    # column: "one" or "double"
    # square: True or False
    # ratio: height/width

    inches_per_pt = 1.0/72.00              # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0     # Most aesthetic ratio
    if (ratio == None): ratio = golden_mean
    if (column == "one"):
        fig_width_pt = columnwidth
    elif (column == "double"):
        fig_width_pt = textwidth
    else:
        raise ValueError("column should be one of ``one'' or ``double''. ")
    fig_width = fig_width_pt*inches_per_pt # Figure width in inches
    if square:
        fig_height = fig_width
    else:
        fig_height = fig_width*ratio
    return [fig_width,fig_height]

errstyle = {'capsize':2, 'ecolor':darkgrey, 'elinewidth':1, 'capthick':1, 'linestyle':'None'}

