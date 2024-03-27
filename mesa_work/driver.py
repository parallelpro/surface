from __future__ import print_function
import re
import numpy as np
import os
import sys
from astropy.io import ascii
import glob
from time import sleep

def set_gyre_inlist(inputFileName, summaryFileName, nu_max, delta_nu):
    # reads in the template inlist and writes out a new inlist with the 
    # parameters set appropriately
    try:
        inlist = open('gyre_template.in','r')
        outlist = open('gyre.in','w')
    except:
        sleep(2.0)
        inlist = open('gyre_template.in','r')
        sleep(2.0)
        outlist = open('gyre.in','w')

    k, b = 0.9638, -1.7145
    width = np.exp(k*np.log(nu_max) + b)
    freqMin = nu_max - width*5 
    freqMax = nu_max + width*5

    # freqMin = nu_max - 10.*delta_nu
    # freqMax = nu_max + 10.*delta_nu
    if freqMin <= 0.: freqMin = 0.0001

    for line in inlist.read().split('\n'):

        first = line.split()
        if len(first)>0:
            if first[0] != '!':

                if re.search("file = 'spb.mesa'",line):
                    line = "\tfile = '%s'  !set by driver.py" % inputFileName

                if re.search("summary\_file = 'summary.txt'",line):
                    line = "\tsummary_file = '%s'  !set by driver.py" % summaryFileName

                if re.search("freq\_min = 100",line):
                    line = "\tfreq_min={:0.5f}  !set by driver.py".format(freqMin)

                if re.search("freq\_max = 1800",line):
                    line = "\tfreq_max={:0.5f}  !set by driver.py".format(freqMax)

        print(line,file=outlist)

    outlist.close()
    inlist.close()


class readTable:
    '''
    A parent class to be wrapped by other class, in order to read in such as mesa history file.
    These files have very similar structures, typically header+table.

    '''
    def __init__(self, filepath, verbose=True):
        '''
        Args:
            filepath: the path of .history file.
            verbose: whether print info to device. default: True.

        Attributes:
            header: a dictionary containing history file header
            track: a structure numpy array containing the evol track
            colnames: a tuple containing column names

        '''  
        self.filepath = filepath
        if verbose: 
            print('Processing :', self.filepath)
        return
    
    def readFile(self, filepath, headerNameLine=1, headerDataLine=2, tableHeaderLine=6):
        '''
        Reads in a file.
        '''

        with open(self.filepath) as f:
            content = [line.split() for line in f]
        header = {content[headerNameLine-1][i]:content[headerDataLine-1][i] for i in range(len(content[headerNameLine-1]))}
        table = np.genfromtxt(self.filepath, skip_header=tableHeaderLine-1, names=True)
        colnames = table.dtype.names

        return header, table, colnames



class history(readTable):
    '''

    A class to read mesa history files, store the data within, and offer useful routines (?).

    '''
    
    def __init__(self, filepath, verbose=True, ifReadProfileIndex=False):
        '''
        Args:
            filepath: the path of .history file.
            verbose: whether print info to device. default: True.

        Attributes:
            header: a dictionary containing history file header
            track: a structure numpy array containing the evol track
            colnames: a tuple containing column names
            profileIndex: a structured array containing the map between model_number and profile_number

        '''
        super().__init__(filepath, verbose)
        self.header, self.track, self.colnames = self.readFile(filepath, headerNameLine=2, headerDataLine=3, tableHeaderLine=6)
        
        if ifReadProfileIndex:
            self.profileIndex = self.read_profile_index()
        return
    
    def read_profile_index(self):
        '''
        Reads in the profile.index file
        '''
        filepath = self.filepath.split('.history')[0] + 'profile.index'
        profileIndex = np.genfromtxt(filepath, skip_header=1, names=('model_number', 'priority', 'profile_number'))
        return profileIndex


class profile(readTable):
    '''

    A class to read mesa history files, store the data within, and offer useful routines.

    '''

    
    def __init__(self, filepath, verbose=True):
        '''
        Args:
            filepath: the path of *profile*.data file.
            verbose: whether print info to device. default: True.

        Attributes:
            header: a dictionary containing history file header
            profile: a structure numpy array containing the structure profile
            colnames: a tuple containing column names

        '''
        super().__init__(filepath, verbose)
        self.header, self.profile, self.colnames = self.readFile(filepath, headerNameLine=2, headerDataLine=3, tableHeaderLine=6)

        return


class sums(readTable):
    '''

    A class to read gyre mode summary file, store the data within, and offer useful routines.

    '''

    
    def __init__(self, filepath, verbose=True):
        '''
        Args:
            filepath: the path of .sums file.
            verbose: whether print info to device. default: True.

        Attributes:
            header: a dictionary containing history file header
            modeSummary: a structure numpy array containing the summary table
            colnames: a tuple containing column names

        '''
        super().__init__(filepath, verbose)
        self.header, self.modeSummary, self.colnames = self.readFile(filepath, headerNameLine=3, headerDataLine=4, tableHeaderLine=6)

        return

def set_mesa_inlist(index, inlist_path, outlist_path, 
                    mass, Xinit, Yinit, Zinit, amlt, 
                    fov_shell, fov0_shell, 
                    fov_core, fov0_core, 
                    ifsetfinalmodel):
    # reads in the template inlist and writes out a new inlist with the 
    # parameters set appropriately
    
    try:
        inlist = open(inlist_path,'r')
        outlist = open(outlist_path,'w')
    except:
        sleep(2.0)
        inlist = open(inlist_path,'r')
        sleep(2.0)
        outlist = open(outlist_path,'w')


    for line in inlist.read().split('\n'):

        first = line.split()
        if len(first)>0:
            if first[0] != '!':

                if re.search('initial\_mass',line):
                    line = "\tinitial_mass = %g  !set by driver.py" % mass

                if re.search('initial\_z =',line):
                    line = "\tinitial_z = %g  !set by driver.py" % Zinit
                
                if re.search('Zbase =',line):
                    line = "\tZbase =%g  !set by driver.py" % Zinit

                if re.search('initial\_y',line):
                    line = "\tinitial_y = %g  !set by driver.py" % Yinit

                if re.search('mixing\_length\_alpha',line):
                    line = "\tmixing_length_alpha = %g  !set by driver.py" % amlt

                if fov_shell >0:
                    if re.search('overshoot\_scheme\(1\)',line):
                        line = "\tovershoot_scheme(1) = '%s' !set by driver.py " % "exponential" 
                    if re.search('overshoot\_zone\_type\(1\)',line):
                        line = "\tovershoot_zone_type(1) = '%s' !set by driver.py " % "any" 
                    if re.search('overshoot\_zone\_loc\(1\)',line):
                        line = "\tovershoot_zone_loc(1) = '%s' !set by driver.py " % "shell" 
                    if re.search('overshoot\_bdy\_loc\(1\)',line):
                        line = "\tovershoot_bdy_loc(1) = '%s' !set by driver.py " % "any" 
                    if re.search('overshoot\_f\(1\)',line):
                        line = "\tovershoot_f(1) = %g  !set by driver.py" %fov_shell
                    if re.search('overshoot\_f0\(1\)',line):
                        line = "\tovershoot_f0(1) = %g !set by driver.py" %fov0_shell

                else:
                    if re.search('overshoot\_scheme\(1\)',line):
                        line = "\t!overshoot_scheme(1) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_zone\_type\(1\)',line):
                        line = "\t!overshoot_zone_type(1) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_zone\_loc\(1\)',line):
                        line = "\t!overshoot_zone_loc(1) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_bdy\_loc\(1\)',line):
                        line = "\t!overshoot_bdy_loc(1) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_f\(1\)',line):
                        line = "\t!overshoot_f(1) = %g  !set by driver.py" % 0.
                    if re.search('overshoot\_f0\(1\)',line):
                        line = "\t!overshoot_f0(1) = %g !set by driver.py" % 0.               

                if fov_core >0:
                    if re.search('overshoot\_scheme\(2\)',line):
                        line = "\tovershoot_scheme(2) = '%s' !set by driver.py " % "exponential" 
                    if re.search('overshoot\_zone\_type\(2\)',line):
                        line = "\tovershoot_zone_type(2) = '%s' !set by driver.py " % "any" 
                    if re.search('overshoot\_zone\_loc\(2\)',line):
                        line = "\tovershoot_zone_loc(2) = '%s' !set by driver.py " % "core" 
                    if re.search('overshoot\_bdy\_loc\(2\)',line):
                        line = "\tovershoot_bdy_loc(2) = '%s' !set by driver.py " % "any" 
                    if re.search('overshoot\_f\(2\)',line):
                        line = "\tovershoot_f(2) = %g  !set by driver.py" %fov_core
                    if re.search('overshoot\_f0\(2\)',line):
                        line = "\tovershoot_f0(2) = %g !set by driver.py" %fov0_core

                else:
                    if re.search('overshoot\_scheme\(2\)',line):
                        line = "\t!overshoot_scheme(2) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_zone\_type\(2\)',line):
                        line = "\t!overshoot_zone_type(2) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_zone\_loc\(2\)',line):
                        line = "\t!overshoot_zone_loc(2) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_bdy\_loc\(2\)',line):
                        line = "\t!overshoot_bdy_loc(2) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_f\(2\)',line):
                        line = "\t!overshoot_f(2) = %g  !set by driver.py" % 0.
                    if re.search('overshoot\_f0\(2\)',line):
                        line = "\t!overshoot_f0(2) = %g !set by driver.py" % 0.   

                # smass = 'm{:03.0f}'.format(mass*100)
                header = 'index{:06.0f}'.format(index)

                if re.search('star\_history\_name',line):
                    output_history_name = header + '.history'
                    line = "\tstar_history_name = '%s'  !set by driver.py" % output_history_name

                if re.search('profile\_data\_prefix',line):
                    output_profile = header + 'profile'
                    line = "\tprofile_data_prefix = '%s' !set by driver.py " % output_profile 
                
                if re.search('profiles\_index\_name =', line):
                    output_index = header + 'profile.index'
                    line = "\tprofiles_index_name = '%s' !set by driver.py " % output_index 

                if ifsetfinalmodel:
                    if re.search('save\_model\_filename =', line):
                        output_final_model = header + 'final.mod'
                        line = "\tsave_model_filename = '%s' !set by driver.py " % output_final_model 

        print(line,file=outlist)

    outlist.close()
    inlist.close()

# free parameters: mass, feh, alpha, y
'''
masses = np.arange(0.4,1.2,0.1) #1.0,0.4
fehs = np.array([-0.4])
alphas = np.array([1.7])
'''

tracks = ascii.read('grid_input_params.txt', delimiter=',')

# ### remove those already existed
# indexes_existed = np.array([int(f[5:11]) for f in os.listdir('../finalmodels/') if f.endswith('final.mod')])
# tracks = tracks[~np.isin(tracks['index'], indexes_existed)]
# ###

Ntracks = len(tracks)
indexes = np.arange(0, Ntracks)
ipart, Npart = int(sys.argv[1]), int(sys.argv[2])

bounds = np.round(np.arange(0,Npart+1)/(Npart)*Ntracks)

idx = (indexes >= bounds[ipart-1]) & (indexes <= bounds[ipart])


for track in tracks[idx]:
    index, mass = track['index'], track['star_mass']
    Xinit, Yinit, Zinit, amlt = track['Xinit'], track['Yinit'], track['Zinit'], track['amlt']
    fov_shell, fov0_shell, fov_core, fov0_core = track['fov_shell'], track['fov0_shell'], track['fov_core'], track['fov0_core']
    mh = np.log10(Zinit/Xinit) - np.log10(0.0134/0.7381) # asplund 09 current abundance
    
    output_history_name = 'index{:06.0f}'.format(index)+'.history'
    output_final_model_name = 'index{:06.0f}final.mod'.format(index)
    output_profile_index_name = 'index{:06.0f}profile.index'.format(index)
    print('Now calculating ',output_history_name)
    # if os.path.exists('LOGS/'+output_history_name): continue
    if os.path.exists('../pre_rgb_tip/finalmodels/'+output_final_model_name): continue

    set_mesa_inlist(index, 'inlist_template', 'inlist', 
                    mass, Xinit, Yinit, Zinit, amlt, 
                    fov_shell, fov0_shell, fov_core, fov0_core, True)
    #os.system('\\rm -r LOGS; \\rm -r png; \\rm -r photos')

    print('------ MESA start ------')
    os.system('./rn1')
    print('------ MESA done ------')

    if os.path.exists(output_final_model_name):
        os.system('mv {:s} ../pre_rgb_tip/finalmodels/'.format(output_final_model_name))

#         # List all files in the home directory
#         files = glob.glob(os.path.expanduser("photos/x*"))
#         if len(files) != 0:
#             # Sort by modification time (mtime) descending
#             latest_binary_photo = sorted(files, key=lambda t: -os.stat(t).st_mtime)[0]
#             os.system('mv {:s} ../finalmodels/{:s}.x'.format(latest_binary_photo, output_final_model_name))
    else:
        os.system('touch ../pre_rgb_tip/finalmodels/{:s}'.format(output_final_model_name))


    # check if history and index file exists
    if (not os.path.exists('LOGS/'+output_history_name)): continue
    if (not os.path.exists('LOGS/'+output_profile_index_name)): continue

    # calculate radial modes
    filepath = 'LOGS/'

    h = history(filepath+output_history_name, ifReadProfileIndex=True)
    track = h.track # 'delta_nu', 'nu_max'
    profileIndex = h.profileIndex # 'model_number', 'priority', 'profile_number'

    fgongPaths = [f for f in os.listdir(filepath) if (f.endswith('.FGONG') & f.startswith('index{:06.0f}'.format(index)))]

    for fgongPath in fgongPaths:

        inputFileName = filepath+fgongPath
        summaryFileName = filepath+fgongPath+'.sum'

        # if os.path.exists(summaryFileName): continue

        profileNo = int(fgongPath.split('profile')[-1].split('.data')[0])
        modelNo = profileIndex['model_number'][profileIndex['profile_number'] == profileNo][0]
        delta_nu = track['delta_nu'][track['model_number'] == modelNo][0]
        nu_max = track['nu_max'][track['model_number'] == modelNo][0]
        
        # correct nu_max with the metallicity dependence
        p = np.array([-0.00718002, -0.05279003,  0.9892853 ])
        nu_max = nu_max * np.polyval(p, mh)
        
        set_gyre_inlist(inputFileName, summaryFileName, nu_max, delta_nu)

        # print('------ GYRE start ------')
        os.system('$GYRE_DIR/bin/gyre gyre.in')
        # print('------ GYRE done ------')


    os.system('mv LOGS/{:s} ../pre_rgb_tip/history/'.format(output_history_name))
    os.system('mv LOGS/{:s} ../pre_rgb_tip/history/'.format(output_profile_index_name))
    if (not os.path.exists('../pre_rgb_tip/freqs/index{:06.0f}/'.format(index))): os.mkdir('../pre_rgb_tip/freqs/index{:06.0f}/'.format(index))
    os.system('mv LOGS/*.sum ../pre_rgb_tip/freqs/index{:06.0f}/'.format(index))


    # os.system('mv LOGS/*.FGONG ../freqs/index{:06.0f}/'.format(index))
    os.system('rm LOGS/*.FGONG')
    os.system('rm LOGS/*profile*.data')


