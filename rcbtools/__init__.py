import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd

__version__ = 0.1
__author__ = "Brad Munson"
__credits__ = "Louisiana State University"

from pkg_resources import resource_filename as rfname
obs_data = rfname(__name__, "observed_abund.csv")

mus = {'H'  : 1.0079,
       'He' : 4.0026,
       'Be' : 9.0122,
       'B'  : 10.811,
       'Li' : 6.941,
       'N'  : 14.0067,
       'O'  : 15.9994,
       'C'  : 12.0107,
       'F'  : 18.9984,
       'Na' : 22.9897,
       'Ne' : 20.1797,
       'Mg' : 24.305,
       'Al' : 26.9815,
       'Si' : 28.0855,
       'P'  : 30.9738,
       'S'  : 32.065,
       'Cl' : 35.453,
       'Ar' : 39.948,
       'K'  : 39.0983,
       'Ca' : 40.078,
       'Ti' : 47.867,
       'V'  : 50.94,
       'Cr' : 51.9961,
       'Mn' : 54.94,
       'Fe' : 55.845,
       'Co' : 58.9332,
       'Ni' : 58.6934,
       'Cu' : 63.546,
       'Zn' : 65.39,
       'Sc' : 44.955908,
       'Y'  : 88.90594,
       'Zr' : 91.224,
       'Ba' : 137.327,
       'La' : 138.90547
       }

def profile2dict(profile, skip_header = 5, global_headers = False):
    '''
    Turn a general profile/history data file into a python dictionary. In a MESA output file,
    the first 5 lines are skipped (default).

    Parameters
    ----------
    profile : str
        A string containing the path to the profile/history data file.
    skip_header : int, optional
        How many lines to skip at beginning of file. The default is 5.
    global_headers : bool, optional
        If True, read the second and third row to get global variables of star.

    Returns
    -------
    p : dict
        Python dictionary containing keys (headers) for data contained in profile/history 
        data file.

    '''
    headers = np.genfromtxt(profile,dtype='str',skip_header=skip_header)[0]
    data = np.loadtxt(profile,skiprows=skip_header+1)
    
    p = dict(zip(headers,data.T))
    if global_headers:
        gheaders = np.genfromtxt(profile,dtype='str',skip_header=1,max_rows = 1)
        gdata = np.loadtxt(profile,skiprows=2,dtype='str',max_rows = 1)
        for i,header in enumerate(gheaders):
            p[header] = gdata[i]
    return p

def makeabund(profile, skip_header = 5, combine_isos = True):
    '''
    Generate a Python Dictionary containing the isotope information from a data file.
    The format of the names must be all lowercase and no non-letter/numeric characters.
    e.g. Hydrogen-1 is h1 and Carbon-12 is c12.
    
    Also generates keys for total abundance of element (adds all isotopes in network).
    e.g. abund['C'] = abund['c8'] + abund['c9'] + ... 

    Parameters
    ----------
    profile : str
        A string containing the path to the profile/history data file.
    return_mus : bool, optional
        If 'True', also return the mean molecular mass dictionary. The default is False.
    skip_header : int, optional
        How many lines to skip at the beginning of file. The default is 5.

    Returns
    -------
    abund : dict
        Dictionary containing mass fraction profile of isotope keys.
    mus : dict, optional
        Dictionary containing mean molecular mass information.

    '''
    headers = np.genfromtxt(profile,dtype='str',skip_header=skip_header)[0]
    data = np.loadtxt(profile,skiprows=skip_header+1)
       
    isos = []
    abunds = {}
    for i,header in enumerate(headers):
        if re.match(r'[a-z]{1,2}[0-9]{1,3}$',header) or header == 'neut':
            isos.append(header)
            abunds[header] = data[:,i]
            
    isos.sort()
    names = [re.sub('[0-9]', '', i) for i in isos]
    for i,name in enumerate(names):
        names[i] = name.capitalize()
            
    
    
    k = 0
    while k < len(isos) and combine_isos:
        eles_idx = np.where(names[k] == np.array(names))[0]
        current_ele = names[k]
        abunds[current_ele] = 0
        for i in eles_idx:
            abunds[current_ele] += abunds[isos[i]]
        k = eles_idx[-1]+1
    return abunds
    
def historyabund(profile, skip_header = 5, combine_isos = True):
    '''
    Take surface abundances from history file formatted from run_star_extras routine.
    Note that the header name needs to be formatted as "rcb_<isotope>". Will also combine
    isotopes for a given element into one surface abundance (Capitalized).

    Parameters
    ----------
    profile : str
        MESA history file where surface abundances can be found.
    skip_header : int, optional
        How many lines to skip in order to reach the data. The default is 5.
    combine_isos : bool, optional
        If True, return combined isotopes for each element in network. The default is True.

    Returns
    -------
    abunds : dict
        Dictionary containing surface abundances at each timestep in history file.

    '''    
    
    headers = np.genfromtxt(profile,dtype='str',skip_header=skip_header)[0]
    data = np.loadtxt(profile,skiprows=skip_header+1)
       
    isos = []
    abunds = {}
    for i,header in enumerate(headers):
        if 'rcb_' in header:
            isos.append(header[4:])
            abunds[header[4:]] = data[:,i]
            
    isos.sort()
    names = [re.sub('[0-9]', '', i) for i in isos]
    for i,name in enumerate(names):
        names[i] = name.capitalize()
            
    
    
    k = 0
    while k < len(isos) and combine_isos:
        eles_idx = np.where(names[k] == np.array(names))[0]
        current_ele = names[k]
        abunds[current_ele] = 0
        for i in eles_idx:
            abunds[current_ele] += abunds[isos[i]]
        k = eles_idx[-1]+1
        

    return abunds

def rcbsurf(*profiles, elements = [], savefig = None, ind_tau = None, labels = [],\
               observed_datafile = obs_data, fig_handles = False, legend = False):
    '''
    This function plots the surface abundances of a list of elements compared to observed RCB stars.
    As opposed to previous iterations, this function will take history datafiles or profiles for any
    list of elements as long as the elements exist in the list of observed RCB stars.
    
    Note: The history datafile can only be used if the surface values are printed to the history log.
    This means you should have a column called "rcb_<isotope>" in your history.data file.

    Parameters
    ----------
    *profiles : str
        Filepaths to profiles or history files.
    elements : list, optional
        List of elements to plot. Format should be uppercase followed by all lowercase.
        e.g. To plot Carbon, Nitrogen, Oxygen, and Lithium use ['C','N','O','Li'].
        
        If elements is left empty, default to ['Li','C','N','O','F','Ne','Na','Mg',
        'Al','Si','S','Ca','Sc','Ti','Ni','Zn','Y','Zr','Ba','La'].
    savefig : str, optional
        If a string is provided, save the plot under that filename. The default is None.
    ind_tau : int, optional
        If an integer is provided, use the first ind_tau rows as the photosphere (this assumes
        that the first row is the stars surface).
        The default is None.
    labels : list, optional
        List of strings for labels for a plot legend. If left empty, the default will be set
        to the profile filenames.
    observed_datafile : str, optional
        String for the filepath of the csv file for the observed abundances.
        The default is 'observed_abund.csv' and should be included in the rcbtools package.
    fig_handles : bool, optional
        If True, return the figure and subplot handles for user control. If you want to change
        the layout margins, first run the command with f, ax = rcbsurf(..., fig_handles = True).
        Then you can change the margins with f.subplots_adjust(). You will need this if using the legend.
    legend : bool, optional
        If True, use "labels" to make a legend on the figure.

    Returns
    -------
    Plots specified abundances.

    '''
    
    if not elements:
        elements = ['Li','C','N','O',\
                    'F','Ne','Na','Mg',\
                    'Al','Si','S','Ca',\
                    'Sc','Ti','Ni','Zn',\
                    'Y','Zr','Ba','La']
        
    file = pd.read_csv(observed_datafile,delimiter=',')
    
    length = len(elements)
    
    row = int((length)**0.5)
    col = int(np.ceil(length/row))
    if row > col: row, col = col, row
    
    f, axarr = plt.subplots(row, col, figsize = [col*3,row*3],sharex = True,sharey=True)
    
    solar = {}
    solar['Li'] = 3.35
    solar['C'] = 8.46
    solar['N'] = 7.90
    solar['O'] = 8.76
    solar['F'] = 4.53
    solar['Ne'] = 7.95
    solar['Na'] = 6.37
    solar['Mg'] = 7.62
    solar['Al'] = 6.54
    solar['Si'] = 7.61
    solar['S'] = 7.26
    solar['Ca'] = 6.41
    solar['Sc'] = 3.15
    solar['Ti'] = 5.00
    solar['Ni'] = 6.29
    solar['Zn'] = 4.70
    solar['Y'] = 2.28
    solar['Zr'] = 2.67
    solar['Ba'] = 2.25
    solar['La'] = 1.25
    solar['Fe'] = 7.54
    
    for i in range(row):
        for j in range(col):
            if col*i+j >= len(elements):
                break
            if elements[col*i+j] in file.columns:
                valid = np.where(file[elements[col*i+j]] == file[elements[col*i+j]])[0]
                axarr[i, j].scatter(file['Fe'][valid]-solar['Fe'], file[elements[col*i+j]][valid]-solar[elements[col*i+j]],\
                                    color=file['Color'][valid],s=100,marker='*')
                axarr[i, j].set_xlim([-2.5,0.0])
                axarr[i, j].xaxis.set_tick_params(labelsize=11)
                axarr[i, j].yaxis.set_tick_params(labelsize=11)
                axarr[i, j].set_ylim([-2.4,3.6])
                axarr[i, j].text(-2.4, 2.5, elements[col*i+j], fontsize=14, fontweight='bold')
                axarr[i, j].yaxis.offsetText.set_visible(False)
                if i == row - 1:
                    axarr[i, j].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
                if j == 0:
                    axarr[i, j].set_ylabel("[X]",fontsize=10, fontweight='bold')
    f.subplots_adjust(left = 0.05, right = 0.95, top = 0.95, bottom = 0.08,\
                        wspace=0, hspace=0)
    
    for index,profile in enumerate(profiles):
        p = profile2dict(profile,global_headers = True)
        if 'zone' in p.keys() or 'history' not in profile:
            abunds = makeabund(profile)
            dm = p['dm']
            if ind_tau:
                ind = ind_tau
            else:
                tau = p['tau']
                ind = np.argmin(abs(tau-1.0))
            for ele in abunds:
                abunds[ele] = sum(dm[0:ind]*abunds[ele][0:ind])/sum(dm[0:ind])
            RCB_age = p['star_age']
            RCB_model = p['model_number']
        else:
            abunds = historyabund(profile)
            ind = np.argmin(p['log_Teff'])
            for ele in abunds:
                abunds[ele] = abunds[ele][ind]
            RCB_age = p['star_age'][ind]
            RCB_model = p['model_number'][ind]
        
        keys = list(abunds.keys())
        
        values = {}
        if 'Fe' in keys:
            values['Fe'] = np.log10(abunds['Fe'])-np.log10(mus['Fe'])+12.15
        else:
            values['Fe'] = solar['Fe']
        print('----------------------------------------------------------------')
        print('PARAMETERS OF INTEREST')
        print('----------------------------------------------------------------')
        print('C12/C13 ', abunds['c12']/abunds['c13'])
        print('O16/O18 ', round(abunds['o16']/abunds['o18'],2))
        print('C/O     ', round(abunds['C']/abunds['O'],2))
        print('RCB Age ', round(float(RCB_age)), 'years')
        print('Model   ', int(RCB_model))
        
        print('----------------------------------------------------------------')	
        print('SURFACE ABUNDANCES FOR SELECTED SPECIES (Solar Unit Abundances)')
        print('----------------------------------------------------------------')
        for element in elements:
            if element in keys:
                values[element] = np.log10(abunds[element])-np.log10(mus[element])+12.15
            else:
                values[element] = -99
            print(element,round(values[element]-solar[element],3))
          
        if not labels:
            labels.append(profile)
        
        for i in range(row):
            for j in range(col):
                if col*i+j >= len(elements):
                    break
                if i == j == 0:
                    axarr[i, j].scatter(values['Fe']-solar['Fe'],\
                                        values[elements[col*i+j]]-solar[elements[col*i+j]],label = labels[index],s=100,marker='s')
                else:
                    axarr[i, j].scatter(values['Fe']-solar['Fe'],\
                                        values[elements[col*i+j]]-solar[elements[col*i+j]],s=100,marker='s')
    
    if legend: f.legend(loc = 'upper right')
    f.align_ylabels(axarr[:,0])
    f.show()
    if savefig:
        f.savefig(savefig)
    if fig_handles: return f, axarr