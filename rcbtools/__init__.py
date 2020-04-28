import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd

__version__ = 0.1
__author__ = "Brad Munson"
__credits__ = "Louisiana State University"

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

def makeabund(profile,return_mus = False, skip_header = 5, combine_isos = True):
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
    
    #Mean atomic masses
    mus = {}
    mus['H']  = 1.0079
    mus['He'] = 4.0026
    mus['Be'] = 9.0122
    mus['B']  = 10.811
    mus['Li'] = 6.941
    mus['N']  = 14.0067
    mus['O']  = 15.9994
    mus['C']  = 12.0107
    mus['F']  = 18.9984
    mus['Na'] = 22.9897
    mus['Ne'] = 20.1797
    mus['Mg'] = 24.305
    mus['Al'] = 26.9815
    mus['Si'] = 28.0855
    mus['P']  = 30.9738
    mus['S']  = 32.065
    mus['Cl'] = 35.453
    mus['Ar'] = 39.948
    mus['K']  = 39.0983
    mus['Ca'] = 40.078
    mus['Ti'] = 47.867
    mus['V']  = 50.94
    mus['Cr'] = 51.9961
    mus['Mn'] = 54.94
    mus['Fe'] = 55.845
    mus['Co'] = 58.9332
    mus['Ni'] = 58.6934
    mus['Cu'] = 63.546
    mus['Zn'] = 65.39
    mus['Sc'] = 44.955908
    mus['Y']  = 88.90594
    mus['Zr'] = 91.224
    mus['Ba'] = 137.327
    mus['La'] = 138.90547
       
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
    if return_mus: 
        return abunds, mus
    else:
        return abunds
    
def historyabund(profile, skip_header = 5, combine_isos = True, return_mus = True):
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
    #Mean atomic masses
    mus = {}
    mus['H']  = 1.0079
    mus['He'] = 4.0026
    mus['Be'] = 9.0122
    mus['B']  = 10.811
    mus['Li'] = 6.941
    mus['N']  = 14.0067
    mus['O']  = 15.9994
    mus['C']  = 12.0107
    mus['F']  = 18.9984
    mus['Na'] = 22.9897
    mus['Ne'] = 20.1797
    mus['Mg'] = 24.305
    mus['Al'] = 26.9815
    mus['Si'] = 28.0855
    mus['P']  = 30.9738
    mus['S']  = 32.065
    mus['Cl'] = 35.453
    mus['Ar'] = 39.948
    mus['K']  = 39.0983
    mus['Ca'] = 40.078
    mus['Ti'] = 47.867
    mus['V']  = 50.94
    mus['Cr'] = 51.9961
    mus['Mn'] = 54.94
    mus['Fe'] = 55.845
    mus['Co'] = 58.9332
    mus['Ni'] = 58.6934
    mus['Cu'] = 63.546
    mus['Zn'] = 65.39
    mus['Sc'] = 44.955908
    mus['Y']  = 88.90594
    mus['Zr'] = 91.224
    mus['Ba'] = 137.327
    mus['La'] = 138.90547
    
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
        
    if return_mus: 
        return abunds, mus
    else:
        return abunds

def surfabund(*profiles, savefig = None, ind_tau = None, labels = []):
    '''
    Creates abundance plots (originally Amber Lauer's plots) for pre-set list isotopes.
    The isotopes in order are 'Li','C','N','O','F','Ne','Na','Mg','Al','Si','S','Ca',
    'Sc','Ti','Ni','Zn','Y','Zr','Ba','La'. You do not need all of these available, if
    the network does not contain one of these elements, it will default to log(X) = -99.

    Parameters
    ----------
    *profiles : str
        Profiles containing the abundance information to plot.
    savefig : str, optional
        If set to a string, save the file under the filename contained in savefig.
        The default is None.
    ind_tau : int, optional
        If a profile does not have a column for the optical depth, specify the index of the row
        that is the lower part of the photosphere from the surface. The default is None.
    labels : list, optional
        List of strings containing the labels for a legend on the plot. The default is [].

    Returns
    -------
    Creates the abundance plot.

    '''
    
    yFull1 = [-1.04,-0.84,-0.74,-1.44,-1.34,-0.64,-1.24,-0.84,-1.04,-0.74,-0.94,-1.24,-1.24,-1.14]
    yFull2 = [-1.94,-1.74,-2.04,-1.54]
    yFull3 = [-1.44,-2.14,-1.24,-0.84,-1.24,-0.54,-0.64,-0.34,-0.74,-0.74,-1.34,-0.74,-0.24,-0.94,-0.94,-0.44]
    
    Li1 = [-0.55,-0.75,0.15,0.15]
    yLi1 = [-1.04,-1.44,-1.24,-0.94]
    
    C1 = [0.74,0.44,0.54,0.34,0.44,0.74,0.14,0.44,0.44,0.44,0.44,0.34,0.34,0.44]
    yC1 = yFull1
    C2 = [0.34,0.34,0.34,1.14]
    yC2 = yFull2
    C3 = [1.04,0.54,0.44,0.54,0.84,0.84,0.64,0.74,0.94,1.24,0.94,1.04,1.04,0.74,0.54,0.94]
    yC3 = yFull3
    
    N1 = [0.50,0.60,1.00,0.60,0.40,0.60,0.40,0.90,0.90,1.20,0.80,0.80,0.80,0.90]
    yN1 = yFull1
    N2 = [0.10,-0.30,0.10,-0.10]
    yN2 = yFull2
    N3 = [0.40,-0.70,0.70,0.60,0.40,0.70,-0.30,0.70,0.30,0.60,0.00,0.40,0.10,0.10,0.50,0.60]
    yN3 = yFull3
    
    O1 = [0.24,-0.86,-0.36,-0.36,0.04,-1.26,-1.06,-0.66,-1.06,-0.36,0.14,-1.06,-0.56,-0.46]
    yO1 = yFull1
    O2 = [-1.26,-0.06,-0.96,0.14]
    yO2 = yFull2
    O3 = [0.54,0.14,-0.26,-0.36,0.14,-0.16,-0.16,-0.66,0.64,0.94,-1.26,0.44,0.14,-0.86,-0.76,-0.46]
    yO3 = yFull3
    
    F1 = [2.37,2.47,2.47,1.67,2.57,2.07,2.67,2.67]
    yF1 = [-1.04,-0.84,-1.44,-0.64,-1.24,-0.84,-1.24,-1.24]
    F2 = [2.17,1.87,1.97]
    yF2 = [-1.94,-1.74,-2.04]
    F3 = [2.77, 1.97,1.97,1.67,2.67]
    yF3 = [-1.44,-2.14,-1.24,-0.84,-1.24]
    
    Ne = [8.3-7.95,7.9-7.95,8.53-7.95,8.62-7.95]
    yNe = [6.5-7.54,5.6-7.54,6.57-7.54,6.47-7.54]
    
    Na1 = [-0.27,-0.07,0.43,-0.77,-0.57,0.03,-0.37,-0.07,-0.07,-0.07,0.03,-0.27,-0.37,-0.37]
    yNa1 = yFull1
    Na2 = [-0.47,-0.57,-0.67,0.03]
    yNa2 = yFull2
    Na3 = [-6.37,-0.97,0.03,-0.07,0.13]
    yNa3 = [-1.44,-2.14,-1.24,-0.84,-1.24]
    
    Mg1 = [-1.92,-1.32,-0.72]
    yMg1 = [-1.24,-1.04,-1.24]
    Mg2 = [-1.52,-1.02,-1.42]
    yMg2 = [-1.94,-2.04,-1.54]
    Mg3 = [-1.62,-1.62,-0.62,-0.72,-0.72,0.18,-0.02,-0.02,-0.02,0.08,-0.42,-0.02,-0.32,-0.32,-0.42,-0.32]
    yMg3 = yFull3
    
    Al1 = [-0.74,-0.54,-1.34,-0.54,-0.24,-0.34,-0.44,-0.34,-0.24,-0.64,-0.84,-0.64]
    yAl1 = [-1.04,-0.84,-1.44,-0.64,-1.24,-0.84,-1.04,-0.74,-0.94,-1.24,-1.24,-1.14]
    Al2 = [-0.94,-1.14,-1.14,-0.84]
    yAl2 = yFull2
    Al3 = [-0.54,-1.84,0.46,-0.54,-1.14,-0.34,-0.84,-0.04,-0.64,-0.34,-0.84,-0.74,-0.64,-0.94,-0.54,-0.24]
    yAl3 = yFull3
    
    Si1 = [-0.41,-0.31,-0.51,-0.91,-0.71,-0.21,-0.61,-0.41,-0.31,-0.21,-0.51,-0.51,-0.41,-0.51]
    ySi1 = yFull1
    Si2 = [-0.11,-0.31,-0.01,-0.61]
    ySi2 = yFull2
    Si3 = [-1.01,-1.51,-0.41,-1.11,-1.71,-0.61,-0.31,-0.51,0.39,0.09,-0.81,-0.41,0.49,0.19,-0.81,-0.51]
    ySi3 = yFull3
    
    S1 = [-0.46,0.04,-0.46,-0.76,-1.06,-0.26,-0.56,-0.36,-0.36,0.44,-0.46,-0.26,-0.26,-0.46]
    yS1 = yFull1
    S2 = [0.14,-0.56,-0.06,-0.86]
    yS2 = yFull2
    S3 = [-0.76,-1.26,-0.16,-0.36,-0.56,-0.06,-0.36,-0.36,-0.16,-0.06,-0.76,-0.26,-0.16,0.54,-0.66,-0.26]
    yS3 = yFull3
    
    Ca1 = [-1.11,-1.11,-1.01,-1.41,-0.91,-0.81,-1.21,-1.01,-1.11,-0.61,-1.01,-1.31,-1.01,-1.11]
    yCa1 = yFull1
    Ca2 = [-1.11,-1.41,-1.21,-1.31]
    yCa2 = yFull2
    Ca3 = [-2.21,-0.81,-0.91,-0.61,-0.21,-0.41,-0.21,-0.61,-1.21,-0.11]
    yCa3 = [-2.14,-1.24,-0.84,-1.24,-0.64,-0.34,-0.74,-0.74,-1.34,-0.74]
    
    Ti1 = [-0.90,-1.00,-1.30,-1.00,-0.90,-0.80]
    yTi1= [-0.84,-0.74,-1.44,-0.64,-1.24,-1.04]
    Ti2 = [-1.50,-1.70,-0.90]
    yTi2 = [-1.94,-2.04,-1.54]
    Ti3 = [-0.70,-1.80,-1.00,-0.40,-0.30,0.20,-0.20,-0.50,-1.10,-0.70]
    yTi3 = [-1.44,-2.14,-1.24,-0.84,-1.24,-0.54,-0.34,-0.74,-1.34,-0.74]
    
    Ni1 = [-0.79,-0.39,-0.19,-0.89,-0.49,-0.09,-0.39,-0.49,-0.29,-0.09,-0.39,-0.49,-0.69,-0.59]
    yNi1 = yFull1
    Ni2 = [-0.49,-1.09,-1.39,-0.39]
    yNi2 = yFull2
    Ni3 = [-2.29,-1.19,-0.59,-0.69,-0.89,-1.19,-0.69]
    yNi3 = [-2.14,-1.24,-0.54,-0.34,-0.74,-1.34,-0.74]
    
    Zn1 = [-0.20,-1.10,0.10,-0.40,-0.30,-0.30,0.00,-0.30,-0.60,-0.30,-0.40]
    yZn1 = [-0.84,-1.44,-0.64,-1.24,-0.84,-1.04,-0.74,-0.94,-1.24,-1.24,-1.14]
    Zn2 = [-0.60,-0.80,-0.80,-0.30]
    yZn2 = yFull2
    Zn3 = [-1.50,-0.20,-0.30]
    yZn3 = [-2.29,-0.89,-1.19]
    
    Sc = np.array([3.2,2.8,2.9,1.9,3.5]) - 3.15
    ySc = np.array([6.3,5.6,6.0,5.4,6.3]) - 7.54
    
    Y = np.array([1.5,1.9,2.0,1.3,1.5,2.8,1.5,2.6,2.4,3.1,2.0,2.0,2.4,1.9,1.3,2.8,1.1,2.2,1.9,2.1,1.9,1.4,2.9,2.2,3.2]) - 2.28
    yY = np.array([6.5,6.7,6.8,6.1,6.2,6.9,6.3,6.7,6.5,6.8,6.6,6.3,6.3,6.4,5.6,5.8,5.5,6.0,6.1,5.4,6.3,6.3,7.0,7.2,6.8]) - 7.54

    Zr = np.array([1.8,1.8,2.3,2.6,1.8,2.3,2.6,1.4,2.1,3.5,1.0,1.9,2.3,3.1,2.7,3.7,3.5]) - 2.67
    yZr = np.array([6.7,6.3,6.7,6.5,6.8,6.3,5.8,5.5,6.0,6.1,5.4,6.3,6.3,7.0,7.2,6.8,6.8]) - 7.54

    Ba = np.array([1.6,1.4,1.5,0.3,1.0,2.1,1.6,2.6,1.5,2.0,1.5,1.6,1.2,1.5,0.9,1.4,0.3,1.3,2.5,0.6,1.7]) - 2.25
    yBa = np.array([6.5,6.7,6.8,6.1,6.2,6.9,6.3,6.7,6.5,6.8,6.6,6.3,6.3,6.4,5.6,5.8,5.5,6.0,6.1,5.4,6.3]) - 7.54

    La = np.array([1.0,1.5,1.5,2.2,1.3,1.9,1.1,0.8,0.5,0.4]) - 1.25
    yLa = np.array([6.7,6.9,6.3,6.7,6.5,6.8,6.6,6.3,5.6,6.0]) - 7.54
    
    f, axarr = plt.subplots(5, 4, figsize=[16,13])
    
    axarr[0, 0].scatter(yLi1, Li1,c='r',s=100,marker='*')
    axarr[0, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[0, 0].set_xlim([-2.5,0.0])
    axarr[0, 0].xaxis.set_tick_params(labelsize=11)
    axarr[0, 0].yaxis.set_tick_params(labelsize=11)
    axarr[0, 0].set_ylim([-3,3.0])
    axarr[0, 0].yaxis.offsetText.set_visible(False)
    axarr[0, 0].text(-2.3, 1.8, 'Li', fontsize=14, fontweight='bold')
    
    axarr[0, 1].scatter(yC1, C1,c='r',s=100,marker='*')
    axarr[0, 1].scatter(yC2, C2,c='b',s=100,marker='*')
    axarr[0, 1].scatter(yC3, C3,c='g',s=100,marker='*')
    axarr[0, 1].set_xlim([-2.5,0.0])
    axarr[0, 1].xaxis.set_tick_params(labelsize=11)
    axarr[0, 1].yaxis.set_tick_params(labelsize=11)
    axarr[0, 1].set_ylim([-3,3.0])
    axarr[0, 1].yaxis.offsetText.set_visible(False)
    axarr[0, 1].text(-2.3, 1.8, 'C', fontsize=14, fontweight='bold')
    
    axarr[0, 2].scatter(yN1, N1,c='r',s=100,marker='*')
    axarr[0, 2].scatter(yN2, N2,c='b',s=100,marker='*')
    axarr[0, 2].scatter(yN3, N3,c='g',s=100,marker='*')
    axarr[0, 2].set_xlim([-2.5,0.0])
    axarr[0, 2].xaxis.set_tick_params(labelsize=11)
    axarr[0, 2].yaxis.set_tick_params(labelsize=11)
    axarr[0, 2].set_ylim([-3,3.0])
    axarr[0, 2].yaxis.offsetText.set_visible(False)
    axarr[0, 2].text(-2.3, 1.8, 'N', fontsize=14, fontweight='bold')
    
    axarr[0, 3].scatter(yO1, O1,c='r',s=100,marker='*')
    axarr[0, 3].scatter(yO2, O2,c='b',s=100,marker='*')
    axarr[0, 3].scatter(yO3, O3,c='g',s=100,marker='*')
    axarr[0, 3].set_xlim([-2.5,0.0])
    axarr[0, 3].xaxis.set_tick_params(labelsize=11)
    axarr[0, 3].yaxis.set_tick_params(labelsize=11)
    axarr[0, 3].set_ylim([-3,3.0])
    axarr[0, 3].yaxis.offsetText.set_visible(False)
    axarr[0, 3].text(-2.3, 1.8, 'O', fontsize=14, fontweight='bold')
    
    axarr[1, 0].scatter(yF1, F1,c='r',s=100,marker='*')
    axarr[1, 0].scatter(yF2, F2,c='b',s=100,marker='*')
    axarr[1, 0].scatter(yF3, F3,c='g',s=100,marker='*')
    axarr[1, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[1, 0].set_xlim([-2.5,0.0])
    axarr[1, 0].xaxis.set_tick_params(labelsize=11)
    axarr[1, 0].yaxis.set_tick_params(labelsize=11)
    axarr[1, 0].set_ylim([-1,4.0])
    axarr[1, 0].yaxis.offsetText.set_visible(False)
    axarr[1, 0].text(-2.3, 3.0, 'F', fontsize=14, fontweight='bold')
    
    axarr[1, 1].scatter(yNe, Ne, c='r',s=100,marker='*')
    axarr[1, 1].set_xlim([-2.5,0.0])
    axarr[1, 1].xaxis.set_tick_params(labelsize=11)
    axarr[1, 1].yaxis.set_tick_params(labelsize=11)
    axarr[1, 1].set_ylim([-3,3.0])
    axarr[1, 1].yaxis.offsetText.set_visible(False)
    axarr[1, 1].text(-2.3, 1.8, 'Ne', fontsize=14, fontweight='bold')
    
    axarr[1, 2].scatter(yNa1, Na1,c='r',s=100,marker='*')
    axarr[1, 2].scatter(yNa2, Na2,c='b',s=100,marker='*')
    axarr[1, 2].scatter(yNa3, Na3,c='g',s=100,marker='*')
    axarr[1, 2].set_xlim([-2.5,0.0])
    axarr[1, 2].xaxis.set_tick_params(labelsize=11)
    axarr[1, 2].yaxis.set_tick_params(labelsize=11)
    axarr[1, 2].set_ylim([-3,3.0])
    axarr[1, 2].yaxis.offsetText.set_visible(False)
    axarr[1, 2].text(-2.3, 1.8, 'Na', fontsize=14, fontweight='bold')
    
    axarr[1, 3].scatter(yMg1, Mg1,c='r',s=100,marker='*')
    axarr[1, 3].scatter(yMg2, Mg2,c='b',s=100,marker='*')
    axarr[1, 3].scatter(yMg3, Mg3,c='g',s=100,marker='*')
    axarr[1, 3].set_xlim([-2.5,0.0])
    axarr[1, 3].xaxis.set_tick_params(labelsize=11)
    axarr[1, 3].yaxis.set_tick_params(labelsize=11)
    axarr[1, 3].set_ylim([-3,3.0])
    axarr[1, 3].yaxis.offsetText.set_visible(False)
    axarr[1, 3].text(-2.3, 1.8, 'Mg', fontsize=14, fontweight='bold')
    
    axarr[2, 0].scatter(yAl1, Al1,c='r',s=100,marker='*')
    axarr[2, 0].scatter(yAl2, Al2,c='b',s=100,marker='*')
    axarr[2, 0].scatter(yAl3, Al3,c='g',s=100,marker='*')
    axarr[2, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[2, 0].set_xlim([-2.5,0.0])
    axarr[2, 0].xaxis.set_tick_params(labelsize=11)
    axarr[2, 0].yaxis.set_tick_params(labelsize=11)
    axarr[2, 0].set_ylim([-3,3.0])
    axarr[2, 0].yaxis.offsetText.set_visible(False)
    axarr[2, 0].text(-2.3, 1.8, 'Al', fontsize=14, fontweight='bold')
    
    axarr[2, 1].scatter(ySi1, Si1,c='r',s=100,marker='*')
    axarr[2, 1].scatter(ySi2, Si2,c='b',s=100,marker='*')
    axarr[2, 1].scatter(ySi3, Si3,c='g',s=100,marker='*')
    axarr[2, 1].set_xlim([-2.5,0.0])
    axarr[2, 1].xaxis.set_tick_params(labelsize=11)
    axarr[2, 1].yaxis.set_tick_params(labelsize=11)
    axarr[2, 1].set_ylim([-3,3.0])
    axarr[2, 1].yaxis.offsetText.set_visible(False)
    axarr[2, 1].text(-2.3, 1.8, 'Si', fontsize=14, fontweight='bold')
    
    axarr[2, 2].scatter(yS1, S1,c='r',s=100,marker='*')
    axarr[2, 2].scatter(yS2, S2,c='b',s=100,marker='*')
    axarr[2, 2].scatter(yS3, S3,c='g',s=100,marker='*')
    axarr[2, 2].set_xlim([-2.5,0.0])
    axarr[2, 2].xaxis.set_tick_params(labelsize=11)
    axarr[2, 2].yaxis.set_tick_params(labelsize=11)
    axarr[2, 2].set_ylim([-3,3.0])
    axarr[2, 2].yaxis.offsetText.set_visible(False)
    axarr[2, 2].text(-2.3, 1.8, 'S', fontsize=14, fontweight='bold')
    
    axarr[2, 3].scatter(yCa1, Ca1,c='r',s=100,marker='*')
    axarr[2, 3].scatter(yCa2, Ca2,c='b',s=100,marker='*')
    axarr[2, 3].scatter(yCa3, Ca3,c='b',s=100,marker='*')
    axarr[2, 3].set_xlim([-2.5,0.0])
    axarr[2, 3].xaxis.set_tick_params(labelsize=11)
    axarr[2, 3].yaxis.set_tick_params(labelsize=11)
    axarr[2, 3].set_ylim([-3,3.0])
    axarr[2, 3].yaxis.offsetText.set_visible(False)
    axarr[2, 3].text(-2.3, 1.8, 'Ca', fontsize=14, fontweight='bold')
    
    axarr[3, 0].scatter(ySc, Sc,c='r',s=100,marker='*')
    axarr[3, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[3, 0].set_xlim([-2.5,0.0])
    axarr[3, 0].xaxis.set_tick_params(labelsize=11)
    axarr[3, 0].yaxis.set_tick_params(labelsize=11)
    axarr[3, 0].set_ylim([-3,3.0])
    axarr[3, 0].yaxis.offsetText.set_visible(False)
    axarr[3, 0].text(-2.3, 1.8, 'Sc', fontsize=14, fontweight='bold')
    
    axarr[3, 1].scatter(yTi1, Ti1,c='r',s=100,marker='*')
    axarr[3, 1].scatter(yTi2, Ti2,c='b',s=100,marker='*')
    axarr[3, 1].scatter(yTi3, Ti3,c='g',s=100,marker='*')
    axarr[3, 1].set_xlim([-2.5,0.0])
    axarr[3, 1].xaxis.set_tick_params(labelsize=11)
    axarr[3, 1].yaxis.set_tick_params(labelsize=11)
    axarr[3, 1].set_ylim([-3,3.0])
    axarr[3, 1].yaxis.offsetText.set_visible(False)
    axarr[3, 1].text(-2.3, 1.8, 'Ti', fontsize=14, fontweight='bold')
    
    axarr[3, 2].scatter(yNi1, Ni1,c='r',s=100,marker='*')
    axarr[3, 2].scatter(yNi2, Ni2,c='b',s=100,marker='*')
    axarr[3, 2].scatter(yNi3, Ni3,c='g',s=100,marker='*')
    axarr[3, 2].set_xlim([-2.5,0.0])
    axarr[3, 2].xaxis.set_tick_params(labelsize=11)
    axarr[3, 2].yaxis.set_tick_params(labelsize=11)
    axarr[3, 2].set_ylim([-3,3.0])
    axarr[3, 2].yaxis.offsetText.set_visible(False)
    axarr[3, 2].text(-2.3, 1.8, 'Ni', fontsize=14, fontweight='bold')
    
    axarr[3, 3].scatter(yZn1, Zn1,c='g',s=100,marker='*')
    axarr[3, 3].scatter(yZn2, Zn2,c='g',s=100,marker='*')
    axarr[3, 3].scatter(yZn3, Zn3,c='g',s=100,marker='*')
    axarr[3, 3].set_xlim([-2.5,0.0])
    axarr[3, 3].xaxis.set_tick_params(labelsize=11)
    axarr[3, 3].yaxis.set_tick_params(labelsize=11)
    axarr[3, 3].set_ylim([-3,3.0])
    axarr[3, 3].yaxis.offsetText.set_visible(False)
    axarr[3, 3].text(-2.3, 1.8, 'Zn', fontsize=14, fontweight='bold')
    
    axarr[4, 0].scatter(yY, Y,c='g',s=100,marker='*')
    axarr[4, 0].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
    axarr[4, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[4, 0].set_xlim([-2.5,0.0])
    axarr[4, 0].xaxis.set_tick_params(labelsize=11)
    axarr[4, 0].yaxis.set_tick_params(labelsize=11)
    axarr[4, 0].set_ylim([-3,3.0])
    axarr[4, 0].yaxis.offsetText.set_visible(False)
    axarr[4, 0].text(-2.3, 1.8, 'Y', fontsize=14, fontweight='bold')
    
    axarr[4, 1].scatter(yZr, Zr,c='g',s=100,marker='*')
    axarr[4, 1].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
    axarr[4, 1].set_xlim([-2.5,0.0])
    axarr[4, 1].xaxis.set_tick_params(labelsize=11)
    axarr[4, 1].yaxis.set_tick_params(labelsize=11)
    axarr[4, 1].set_ylim([-3,3.0])
    axarr[4, 1].yaxis.offsetText.set_visible(False)
    axarr[4, 1].text(-2.3, 1.8, 'Zr', fontsize=14, fontweight='bold')
    
    axarr[4, 2].scatter(yBa, Ba,c='g',s=100,marker='*')
    axarr[4, 2].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
    axarr[4, 2].set_xlim([-2.5,0.0])
    axarr[4, 2].xaxis.set_tick_params(labelsize=11)
    axarr[4, 2].yaxis.set_tick_params(labelsize=11)
    axarr[4, 2].set_ylim([-3,3.0])
    axarr[4, 2].yaxis.offsetText.set_visible(False)
    axarr[4, 2].text(-2.3, 1.8, 'Ba', fontsize=14, fontweight='bold')
    
    axarr[4, 3].scatter(yLa, La,c='g',s=100,marker='*')
    axarr[4, 3].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
    axarr[4, 3].set_xlim([-2.5,0.0])
    axarr[4, 3].xaxis.set_tick_params(labelsize=11)
    axarr[4, 3].yaxis.set_tick_params(labelsize=11)
    axarr[4, 3].set_ylim([-3,3.0])
    axarr[4, 3].yaxis.offsetText.set_visible(False)
    axarr[4, 3].text(-2.3, 1.8, 'La', fontsize=14, fontweight='bold')
    
    
    for index,profile in enumerate(profiles):
        abunds,mus = makeabund(profile,return_mus=True)
        p = profile2dict(profile)
        
        dm = p['dm']
        if ind_tau:
            ind = ind_tau
        else:
            tau = p['tau']
            ind = np.argmin(abs(tau-1.0))
            
        o16tot = sum(dm[0:ind]*abunds['o16'][0:ind])/sum(dm[0:ind])
        o18tot = sum(dm[0:ind]*abunds['o18'][0:ind])/sum(dm[0:ind])
        c12tot = sum(dm[0:ind]*abunds['c12'][0:ind])/sum(dm[0:ind])
        c13tot = sum(dm[0:ind]*abunds['c13'][0:ind])/sum(dm[0:ind])
        
        keys = list(abunds.keys())
        tots = list(mus.keys())
        
        values = np.zeros(len(tots))
        print('----------------------------------------------------------------')	
        print('SURFACE ABUNDANCES FOR SELECTED SPECIES (TOTAL FOR EACH ELEMENT)')
        print('----------------------------------------------------------------')
        for i,tot in enumerate(tots):
            if tot in keys:
                values[i] = np.log10(sum(dm[0:ind]*abunds[tot][0:ind])/sum(dm[0:ind]))-np.log10(mus[tot])+12.15
                print(tot, ' ', round(values[i],2))
        print('----------------------------------------------------------------')
        print('SURFACE ABUNDANCE RATIOS FOR ELEMENTS OF INTEREST')
        print('----------------------------------------------------------------')
        print('C12/C13 ', c12tot/c13tot)
        print('O16/O18 ', round(o16tot/o18tot,2))
        print('C/O     ', round(sum(dm[0:ind]*abunds['C'][0:ind]/sum(dm[0:ind])/(sum(dm[0:ind]*abunds['O'][0:ind])/sum(dm[0:ind])),2)))
        
        tots_name = np.array(['Li','C','N','O',\
                              'F','Ne','Na','Mg',\
                              'Al','Si','S','Ca',\
                              'Sc','Ti','Ni','Zn',\
                              'Y','Zr','Ba','La',\
                                  'Fe'])
        tots_val = np.zeros(len(tots_name))-99
        tots_val[12] = 7.54
        
        for j,tot in enumerate(tots):
            for i,name in enumerate(tots_name):
                if name == tot and tot in keys:
                    tots_val[i] = values[j]
        
        Lodders = np.array([3.35,8.46,7.90,8.76,\
                            4.53,7.95,6.37,7.62,\
                            6.54,7.61,7.26,6.41,\
                            3.15,5.00,6.29,4.70,\
                            2.28,2.67,2.25,1.25,\
                                7.54])
        ele = tots_val - Lodders
        
        print('----------------------------------------------------------------')	
        print('SURFACE ABUNDANCES FOR SELECTED SPECIES (Solar Unit Abundances)')
        print('----------------------------------------------------------------')
        for idx,name in enumerate(tots_name):
            print(name,round(ele[idx],2))
            
        if not labels:
            labels.append(profile)
        
        axarr[0, 0].scatter(ele[20],ele[0],label=labels[index],s=100,marker='s')
        axarr[0, 1].scatter(ele[20],ele[1],s=100,marker='s')
        axarr[0, 2].scatter(ele[20],ele[2],s=100,marker='s')
        axarr[0, 3].scatter(ele[20],ele[3],s=100,marker='s')
        axarr[1, 0].scatter(ele[20],ele[4],s=100,marker='s')
        axarr[1, 1].scatter(ele[20],ele[5],s=100,marker='s')
        axarr[1, 2].scatter(ele[20],ele[6],s=100,marker='s')
        axarr[1, 3].scatter(ele[20],ele[7],s=100,marker='s')
        axarr[2, 0].scatter(ele[20],ele[8],s=100,marker='s')
        axarr[2, 1].scatter(ele[20],ele[9],s=100,marker='s')
        axarr[2, 2].scatter(ele[20],ele[10],s=100,marker='s')
        axarr[2, 3].scatter(ele[20],ele[11],s=100,marker='s')
        axarr[3, 0].scatter(ele[20],ele[12],s=100,marker='s')
        axarr[3, 1].scatter(ele[20],ele[13],s=100,marker='s')
        axarr[3, 2].scatter(ele[20],ele[14],s=100,marker='s')
        axarr[3, 3].scatter(ele[20],ele[15],s=100,marker='s')
        axarr[4, 0].scatter(ele[20],ele[16],s=100,marker='s')
        axarr[4, 1].scatter(ele[20],ele[17],s=100,marker='s')
        axarr[4, 2].scatter(ele[20],ele[18],s=100,marker='s')
        axarr[4, 3].scatter(ele[20],ele[19],s=100,marker='s')
    
    f.legend(loc = 'upper left')
    f.align_ylabels(axarr[:,0])
    plt.show()
    if savefig:
        f.savefig(savefig)

def surfabund2(*profiles, elements, savefig = None, ind_tau = None, labels = [],\
               observed_datafile = 'rcbtools.__path__/observed_abund.csv'):
    '''
    Similar to surfabund, but more general. The user will provide a list of elements to plot.
    Also, the star symbols representing observed surface abundances are automatically plotted
    using a csv file specified by the observed_datafile argument. 

    Parameters
    ----------
    *profiles : str
        Strings of filepaths containing isotope profiles.
    elements : list
        List of elements to plot. Format should be uppercase followed by all lowercase.
        e.g. To plot Carbon, Nitrogen, Oxygen, and Lithium use ['C','N','O','Li'].
    savefig : str, optional
        If a string is provided, save the plot under that filename. The default is None.
    ind_tau : int, optional
        If an integer is provided, use the first ind_tau rows as the photosphere.
        The default is None.
    labels : list, optional
        List of strings for labels for a plot legend. The default is [].
    observed_datafile : str, optional
        String for the filepath of the csv file for the observed abundances.
        The default is 'observed_abund.csv'.

    Returns
    -------
    Plots specified abundances.

    '''
    
    file = pd.read_csv(observed_datafile,delimiter=',')
    
    length = len(elements)
    
    row = int((length)**0.5)
    col = int(np.ceil(length/row))
    if row > col: row, col = col, row
    
    f, axarr = plt.subplots(row, col, figsize = [col*3,row*3])
    
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
            if elements[col*i+j] in file.columns:
                axarr[i, j].scatter(file['Fe']-solar['Fe'], file[elements[col*i+j]]-solar[elements[col*i+j]],\
                                    c='r',s=100,marker='*')
                axarr[i, j].set_xlim([-2.5,0.0])
                axarr[i, j].xaxis.set_tick_params(labelsize=11)
                axarr[i, j].yaxis.set_tick_params(labelsize=11)
                if elements[col*i+j] == 'F':
                    axarr[i, j].set_ylim([-1,4.0])
                    axarr[i, j].text(-2.3, 3.0, elements[col*i+j], fontsize=14, fontweight='bold')
                else:
                    axarr[i, j].set_ylim([-3,3.0])
                    axarr[i, j].text(-2.3, 1.8, elements[col*i+j], fontsize=14, fontweight='bold')
                axarr[i, j].yaxis.offsetText.set_visible(False)
                if i == row - 1:
                    axarr[i, j].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
                if j == 0:
                    axarr[i, j].set_ylabel("[X]",fontsize=10, fontweight='bold')
    
    
    for index,profile in enumerate(profiles):
        abunds,mus = makeabund(profile,return_mus=True)
        p = profile2dict(profile)
        
        dm = p['dm']
        if ind_tau:
            ind = ind_tau
        else:
            tau = p['tau']
            ind = np.argmin(abs(tau-1.0))
            
        o16tot = sum(dm[0:ind]*abunds['o16'][0:ind])/sum(dm[0:ind])
        o18tot = sum(dm[0:ind]*abunds['o18'][0:ind])/sum(dm[0:ind])
        c12tot = sum(dm[0:ind]*abunds['c12'][0:ind])/sum(dm[0:ind])
        c13tot = sum(dm[0:ind]*abunds['c13'][0:ind])/sum(dm[0:ind])
        
        keys = list(abunds.keys())
        
        values = {}
        if 'Fe' in keys:
            values['Fe'] = np.log10(sum(dm[0:ind]*abunds['Fe'][0:ind])/sum(dm[0:ind]))\
                                -np.log10(mus['Fe'])+12.15
        else:
            values['Fe'] = solar['Fe']
        print('----------------------------------------------------------------')	
        print('SURFACE ABUNDANCES FOR SELECTED SPECIES (TOTAL FOR EACH ELEMENT)')
        print('----------------------------------------------------------------')
        for element in elements:
            if element in keys:
                values[element] = np.log10(sum(dm[0:ind]*abunds[element][0:ind])/sum(dm[0:ind]))\
                                -np.log10(mus[element])+12.15
            else:
                values[element] = -99
            print(element, ' ', round(values[element],3))
            
        print('----------------------------------------------------------------')
        print('SURFACE ABUNDANCE RATIOS FOR ELEMENTS OF INTEREST')
        print('----------------------------------------------------------------')
        print('C12/C13 ', c12tot/c13tot)
        print('O16/O18 ', round(o16tot/o18tot,3))
        print('C/O     ', round(sum(dm[0:ind]*abunds['C'][0:ind]/sum(dm[0:ind])/(sum(dm[0:ind]*abunds['O'][0:ind])/sum(dm[0:ind])),3)))
        
        print('----------------------------------------------------------------')	
        print('SURFACE ABUNDANCES FOR SELECTED SPECIES (Solar Unit Abundances)')
        print('----------------------------------------------------------------')
        for element in elements:
            print(element,round(values[element]-solar[element],3))
          
        if not labels:
            labels.append(profile)
        
        for i in range(row):
            for j in range(col):
                if i == j == 0:
                    axarr[i, j].scatter(values['Fe']-solar['Fe'],\
                                        values[elements[col*i+j]]-solar[elements[col*i+j]],label = labels[index],s=100,marker='s')
                else:
                    axarr[i, j].scatter(values['Fe']-solar['Fe'],\
                                        values[elements[col*i+j]]-solar[elements[col*i+j]],s=100,marker='s')
    
    f.legend(loc = 'upper left')
    f.align_ylabels(axarr[:,0])
    plt.show()
    if savefig:
        f.savefig(savefig)
        
def surfabund_hist(*filenames, labels):
    '''
    Generate a plot of key surface abundances and observational values for model.
    Input argument "filenames" should be history.data files from MESA LOGS/ directories.
    Defines RCB phase by model with the minimum log_Teff.
    
    Note: This only works if MESA has history columns "rcb_<isotope>" as an output.

    Parameters
    ----------
    *filenames : str
        history.data files used for analysis
    labels : list
        List of labels used for legend in plot.

    Returns
    -------
    An array of subplots containing information of surface abundance for key species
    at the RCB phase.

    '''
    
    
    yFull1 = [-1.04,-0.84,-0.74,-1.44,-1.34,-0.64,-1.24,-0.84,-1.04,-0.74,-0.94,-1.24,-1.24,-1.14]
    yFull2 = [-1.94,-1.74,-2.04,-1.54]
    yFull3 = [-1.44,-2.14,-1.24,-0.84,-1.24,-0.54,-0.64,-0.34,-0.74,-0.74,-1.34,-0.74,-0.24,-0.94,-0.94,-0.44]
    
    Li1 = [-0.55,-0.75,0.15,0.15]
    yLi1 = [-1.04,-1.44,-1.24,-0.94]
    
    C1 = [0.74,0.44,0.54,0.34,0.44,0.74,0.14,0.44,0.44,0.44,0.44,0.34,0.34,0.44]
    yC1 = yFull1
    C2 = [0.34,0.34,0.34,1.14]
    yC2 = yFull2
    C3 = [1.04,0.54,0.44,0.54,0.84,0.84,0.64,0.74,0.94,1.24,0.94,1.04,1.04,0.74,0.54,0.94]
    yC3 = yFull3
    
    N1 = [0.50,0.60,1.00,0.60,0.40,0.60,0.40,0.90,0.90,1.20,0.80,0.80,0.80,0.90]
    yN1 = yFull1
    N2 = [0.10,-0.30,0.10,-0.10]
    yN2 = yFull2
    N3 = [0.40,-0.70,0.70,0.60,0.40,0.70,-0.30,0.70,0.30,0.60,0.00,0.40,0.10,0.10,0.50,0.60]
    yN3 = yFull3
    
    O1 = [0.24,-0.86,-0.36,-0.36,0.04,-1.26,-1.06,-0.66,-1.06,-0.36,0.14,-1.06,-0.56,-0.46]
    yO1 = yFull1
    O2 = [-1.26,-0.06,-0.96,0.14]
    yO2 = yFull2
    O3 = [0.54,0.14,-0.26,-0.36,0.14,-0.16,-0.16,-0.66,0.64,0.94,-1.26,0.44,0.14,-0.86,-0.76,-0.46]
    yO3 = yFull3
    
    F1 = [2.37,2.47,2.47,1.67,2.57,2.07,2.67,2.67]
    yF1 = [-1.04,-0.84,-1.44,-0.64,-1.24,-0.84,-1.24,-1.24]
    F2 = [2.17,1.87,1.97]
    yF2 = [-1.94,-1.74,-2.04]
    F3 = [2.77, 1.97,1.97,1.67,2.67]
    yF3 = [-1.44,-2.14,-1.24,-0.84,-1.24]
    
    Ne = [8.3-7.95,7.9-7.95,8.53-7.95,8.62-7.95]
    yNe = [6.5-7.54,5.6-7.54,6.57-7.54,6.47-7.54]
    
    Na1 = [-0.27,-0.07,0.43,-0.77,-0.57,0.03,-0.37,-0.07,-0.07,-0.07,0.03,-0.27,-0.37,-0.37]
    yNa1 = yFull1
    Na2 = [-0.47,-0.57,-0.67,0.03]
    yNa2 = yFull2
    Na3 = [-6.37,-0.97,0.03,-0.07,0.13]
    yNa3 = [-1.44,-2.14,-1.24,-0.84,-1.24]
    
    Mg1 = [-1.92,-1.32,-0.72]
    yMg1 = [-1.24,-1.04,-1.24]
    Mg2 = [-1.52,-1.02,-1.42]
    yMg2 = [-1.94,-2.04,-1.54]
    Mg3 = [-1.62,-1.62,-0.62,-0.72,-0.72,0.18,-0.02,-0.02,-0.02,0.08,-0.42,-0.02,-0.32,-0.32,-0.42,-0.32]
    yMg3 = yFull3
    
    Al1 = [-0.74,-0.54,-1.34,-0.54,-0.24,-0.34,-0.44,-0.34,-0.24,-0.64,-0.84,-0.64]
    yAl1 = [-1.04,-0.84,-1.44,-0.64,-1.24,-0.84,-1.04,-0.74,-0.94,-1.24,-1.24,-1.14]
    Al2 = [-0.94,-1.14,-1.14,-0.84]
    yAl2 = yFull2
    Al3 = [-0.54,-1.84,0.46,-0.54,-1.14,-0.34,-0.84,-0.04,-0.64,-0.34,-0.84,-0.74,-0.64,-0.94,-0.54,-0.24]
    yAl3 = yFull3
    
    Si1 = [-0.41,-0.31,-0.51,-0.91,-0.71,-0.21,-0.61,-0.41,-0.31,-0.21,-0.51,-0.51,-0.41,-0.51]
    ySi1 = yFull1
    Si2 = [-0.11,-0.31,-0.01,-0.61]
    ySi2 = yFull2
    Si3 = [-1.01,-1.51,-0.41,-1.11,-1.71,-0.61,-0.31,-0.51,0.39,0.09,-0.81,-0.41,0.49,0.19,-0.81,-0.51]
    ySi3 = yFull3
    
    S1 = [-0.46,0.04,-0.46,-0.76,-1.06,-0.26,-0.56,-0.36,-0.36,0.44,-0.46,-0.26,-0.26,-0.46]
    yS1 = yFull1
    S2 = [0.14,-0.56,-0.06,-0.86]
    yS2 = yFull2
    S3 = [-0.76,-1.26,-0.16,-0.36,-0.56,-0.06,-0.36,-0.36,-0.16,-0.06,-0.76,-0.26,-0.16,0.54,-0.66,-0.26]
    yS3 = yFull3
    
    Ca1 = [-1.11,-1.11,-1.01,-1.41,-0.91,-0.81,-1.21,-1.01,-1.11,-0.61,-1.01,-1.31,-1.01,-1.11]
    yCa1 = yFull1
    Ca2 = [-1.11,-1.41,-1.21,-1.31]
    yCa2 = yFull2
    Ca3 = [-2.21,-0.81,-0.91,-0.61,-0.21,-0.41,-0.21,-0.61,-1.21,-0.11]
    yCa3 = [-2.14,-1.24,-0.84,-1.24,-0.64,-0.34,-0.74,-0.74,-1.34,-0.74]
    
    Ti1 = [-0.90,-1.00,-1.30,-1.00,-0.90,-0.80]
    yTi1= [-0.84,-0.74,-1.44,-0.64,-1.24,-1.04]
    Ti2 = [-1.50,-1.70,-0.90]
    yTi2 = [-1.94,-2.04,-1.54]
    Ti3 = [-0.70,-1.80,-1.00,-0.40,-0.30,0.20,-0.20,-0.50,-1.10,-0.70]
    yTi3 = [-1.44,-2.14,-1.24,-0.84,-1.24,-0.54,-0.34,-0.74,-1.34,-0.74]
    
    Ni1 = [-0.79,-0.39,-0.19,-0.89,-0.49,-0.09,-0.39,-0.49,-0.29,-0.09,-0.39,-0.49,-0.69,-0.59]
    yNi1 = yFull1
    Ni2 = [-0.49,-1.09,-1.39,-0.39]
    yNi2 = yFull2
    Ni3 = [-2.29,-1.19,-0.59,-0.69,-0.89,-1.19,-0.69]
    yNi3 = [-2.14,-1.24,-0.54,-0.34,-0.74,-1.34,-0.74]
    
    Zn1 = [-0.20,-1.10,0.10,-0.40,-0.30,-0.30,0.00,-0.30,-0.60,-0.30,-0.40]
    yZn1 = [-0.84,-1.44,-0.64,-1.24,-0.84,-1.04,-0.74,-0.94,-1.24,-1.24,-1.14]
    Zn2 = [-0.60,-0.80,-0.80,-0.30]
    yZn2 = yFull2
    Zn3 = [-1.50,-0.20,-0.30]
    yZn3 = [-2.29,-0.89,-1.19]
    
    Sc = np.array([3.2,2.8,2.9,1.9,3.5]) - 3.15
    ySc = np.array([6.3,5.6,6.0,5.4,6.3]) - 7.54
    
    Y = np.array([1.5,1.9,2.0,1.3,1.5,2.8,1.5,2.6,2.4,3.1,2.0,2.0,2.4,1.9,1.3,2.8,1.1,2.2,1.9,2.1,1.9,1.4,2.9,2.2,3.2]) - 2.28
    yY = np.array([6.5,6.7,6.8,6.1,6.2,6.9,6.3,6.7,6.5,6.8,6.6,6.3,6.3,6.4,5.6,5.8,5.5,6.0,6.1,5.4,6.3,6.3,7.0,7.2,6.8]) - 7.54

    Zr = np.array([1.8,1.8,2.3,2.6,1.8,2.3,2.6,1.4,2.1,3.5,1.0,1.9,2.3,3.1,2.7,3.7,3.5]) - 2.67
    yZr = np.array([6.7,6.3,6.7,6.5,6.8,6.3,5.8,5.5,6.0,6.1,5.4,6.3,6.3,7.0,7.2,6.8,6.8]) - 7.54

    Ba = np.array([1.6,1.4,1.5,0.3,1.0,2.1,1.6,2.6,1.5,2.0,1.5,1.6,1.2,1.5,0.9,1.4,0.3,1.3,2.5,0.6,1.7]) - 2.25
    yBa = np.array([6.5,6.7,6.8,6.1,6.2,6.9,6.3,6.7,6.5,6.8,6.6,6.3,6.3,6.4,5.6,5.8,5.5,6.0,6.1,5.4,6.3]) - 7.54

    La = np.array([1.0,1.5,1.5,2.2,1.3,1.9,1.1,0.8,0.5,0.4]) - 1.25
    yLa = np.array([6.7,6.9,6.3,6.7,6.5,6.8,6.6,6.3,5.6,6.0]) - 7.54
    
    f, axarr = plt.subplots(5, 4, figsize=[16,13])
    
    axarr[0, 0].scatter(yLi1, Li1,c='r',s=100,marker='*')
    axarr[0, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[0, 0].set_xlim([-2.5,0.0])
    axarr[0, 0].xaxis.set_tick_params(labelsize=11)
    axarr[0, 0].yaxis.set_tick_params(labelsize=11)
    axarr[0, 0].set_ylim([-3,3.0])
    axarr[0, 0].yaxis.offsetText.set_visible(False)
    axarr[0, 0].text(-2.3, 1.8, 'Li', fontsize=14, fontweight='bold')
    
    axarr[0, 1].scatter(yC1, C1,c='r',s=100,marker='*')
    axarr[0, 1].scatter(yC2, C2,c='b',s=100,marker='*')
    axarr[0, 1].scatter(yC3, C3,c='g',s=100,marker='*')
    axarr[0, 1].set_xlim([-2.5,0.0])
    axarr[0, 1].xaxis.set_tick_params(labelsize=11)
    axarr[0, 1].yaxis.set_tick_params(labelsize=11)
    axarr[0, 1].set_ylim([-3,3.0])
    axarr[0, 1].yaxis.offsetText.set_visible(False)
    axarr[0, 1].text(-2.3, 1.8, 'C', fontsize=14, fontweight='bold')
    
    axarr[0, 2].scatter(yN1, N1,c='r',s=100,marker='*')
    axarr[0, 2].scatter(yN2, N2,c='b',s=100,marker='*')
    axarr[0, 2].scatter(yN3, N3,c='g',s=100,marker='*')
    axarr[0, 2].set_xlim([-2.5,0.0])
    axarr[0, 2].xaxis.set_tick_params(labelsize=11)
    axarr[0, 2].yaxis.set_tick_params(labelsize=11)
    axarr[0, 2].set_ylim([-3,3.0])
    axarr[0, 2].yaxis.offsetText.set_visible(False)
    axarr[0, 2].text(-2.3, 1.8, 'N', fontsize=14, fontweight='bold')
    
    axarr[0, 3].scatter(yO1, O1,c='r',s=100,marker='*')
    axarr[0, 3].scatter(yO2, O2,c='b',s=100,marker='*')
    axarr[0, 3].scatter(yO3, O3,c='g',s=100,marker='*')
    axarr[0, 3].set_xlim([-2.5,0.0])
    axarr[0, 3].xaxis.set_tick_params(labelsize=11)
    axarr[0, 3].yaxis.set_tick_params(labelsize=11)
    axarr[0, 3].set_ylim([-3,3.0])
    axarr[0, 3].yaxis.offsetText.set_visible(False)
    axarr[0, 3].text(-2.3, 1.8, 'O', fontsize=14, fontweight='bold')
    
    axarr[1, 0].scatter(yF1, F1,c='r',s=100,marker='*')
    axarr[1, 0].scatter(yF2, F2,c='b',s=100,marker='*')
    axarr[1, 0].scatter(yF3, F3,c='g',s=100,marker='*')
    axarr[1, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[1, 0].set_xlim([-2.5,0.0])
    axarr[1, 0].xaxis.set_tick_params(labelsize=11)
    axarr[1, 0].yaxis.set_tick_params(labelsize=11)
    axarr[1, 0].set_ylim([-1,4.0])
    axarr[1, 0].yaxis.offsetText.set_visible(False)
    axarr[1, 0].text(-2.3, 3.0, 'F', fontsize=14, fontweight='bold')
    
    axarr[1, 1].scatter(yNe, Ne, c='r',s=100,marker='*')
    axarr[1, 1].set_xlim([-2.5,0.0])
    axarr[1, 1].xaxis.set_tick_params(labelsize=11)
    axarr[1, 1].yaxis.set_tick_params(labelsize=11)
    axarr[1, 1].set_ylim([-3,3.0])
    axarr[1, 1].yaxis.offsetText.set_visible(False)
    axarr[1, 1].text(-2.3, 1.8, 'Ne', fontsize=14, fontweight='bold')
    
    axarr[1, 2].scatter(yNa1, Na1,c='r',s=100,marker='*')
    axarr[1, 2].scatter(yNa2, Na2,c='b',s=100,marker='*')
    axarr[1, 2].scatter(yNa3, Na3,c='g',s=100,marker='*')
    axarr[1, 2].set_xlim([-2.5,0.0])
    axarr[1, 2].xaxis.set_tick_params(labelsize=11)
    axarr[1, 2].yaxis.set_tick_params(labelsize=11)
    axarr[1, 2].set_ylim([-3,3.0])
    axarr[1, 2].yaxis.offsetText.set_visible(False)
    axarr[1, 2].text(-2.3, 1.8, 'Na', fontsize=14, fontweight='bold')
    
    axarr[1, 3].scatter(yMg1, Mg1,c='r',s=100,marker='*')
    axarr[1, 3].scatter(yMg2, Mg2,c='b',s=100,marker='*')
    axarr[1, 3].scatter(yMg3, Mg3,c='g',s=100,marker='*')
    axarr[1, 3].set_xlim([-2.5,0.0])
    axarr[1, 3].xaxis.set_tick_params(labelsize=11)
    axarr[1, 3].yaxis.set_tick_params(labelsize=11)
    axarr[1, 3].set_ylim([-3,3.0])
    axarr[1, 3].yaxis.offsetText.set_visible(False)
    axarr[1, 3].text(-2.3, 1.8, 'Mg', fontsize=14, fontweight='bold')
    
    axarr[2, 0].scatter(yAl1, Al1,c='r',s=100,marker='*')
    axarr[2, 0].scatter(yAl2, Al2,c='b',s=100,marker='*')
    axarr[2, 0].scatter(yAl3, Al3,c='g',s=100,marker='*')
    axarr[2, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[2, 0].set_xlim([-2.5,0.0])
    axarr[2, 0].xaxis.set_tick_params(labelsize=11)
    axarr[2, 0].yaxis.set_tick_params(labelsize=11)
    axarr[2, 0].set_ylim([-3,3.0])
    axarr[2, 0].yaxis.offsetText.set_visible(False)
    axarr[2, 0].text(-2.3, 1.8, 'Al', fontsize=14, fontweight='bold')
    
    axarr[2, 1].scatter(ySi1, Si1,c='r',s=100,marker='*')
    axarr[2, 1].scatter(ySi2, Si2,c='b',s=100,marker='*')
    axarr[2, 1].scatter(ySi3, Si3,c='g',s=100,marker='*')
    axarr[2, 1].set_xlim([-2.5,0.0])
    axarr[2, 1].xaxis.set_tick_params(labelsize=11)
    axarr[2, 1].yaxis.set_tick_params(labelsize=11)
    axarr[2, 1].set_ylim([-3,3.0])
    axarr[2, 1].yaxis.offsetText.set_visible(False)
    axarr[2, 1].text(-2.3, 1.8, 'Si', fontsize=14, fontweight='bold')
    
    axarr[2, 2].scatter(yS1, S1,c='r',s=100,marker='*')
    axarr[2, 2].scatter(yS2, S2,c='b',s=100,marker='*')
    axarr[2, 2].scatter(yS3, S3,c='g',s=100,marker='*')
    axarr[2, 2].set_xlim([-2.5,0.0])
    axarr[2, 2].xaxis.set_tick_params(labelsize=11)
    axarr[2, 2].yaxis.set_tick_params(labelsize=11)
    axarr[2, 2].set_ylim([-3,3.0])
    axarr[2, 2].yaxis.offsetText.set_visible(False)
    axarr[2, 2].text(-2.3, 1.8, 'S', fontsize=14, fontweight='bold')
    
    axarr[2, 3].scatter(yCa1, Ca1,c='r',s=100,marker='*')
    axarr[2, 3].scatter(yCa2, Ca2,c='b',s=100,marker='*')
    axarr[2, 3].scatter(yCa3, Ca3,c='b',s=100,marker='*')
    axarr[2, 3].set_xlim([-2.5,0.0])
    axarr[2, 3].xaxis.set_tick_params(labelsize=11)
    axarr[2, 3].yaxis.set_tick_params(labelsize=11)
    axarr[2, 3].set_ylim([-3,3.0])
    axarr[2, 3].yaxis.offsetText.set_visible(False)
    axarr[2, 3].text(-2.3, 1.8, 'Ca', fontsize=14, fontweight='bold')
    
    axarr[3, 0].scatter(ySc, Sc,c='r',s=100,marker='*')
    axarr[3, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[3, 0].set_xlim([-2.5,0.0])
    axarr[3, 0].xaxis.set_tick_params(labelsize=11)
    axarr[3, 0].yaxis.set_tick_params(labelsize=11)
    axarr[3, 0].set_ylim([-3,3.0])
    axarr[3, 0].yaxis.offsetText.set_visible(False)
    axarr[3, 0].text(-2.3, 1.8, 'Sc', fontsize=14, fontweight='bold')
    
    axarr[3, 1].scatter(yTi1, Ti1,c='r',s=100,marker='*')
    axarr[3, 1].scatter(yTi2, Ti2,c='b',s=100,marker='*')
    axarr[3, 1].scatter(yTi3, Ti3,c='g',s=100,marker='*')
    axarr[3, 1].set_xlim([-2.5,0.0])
    axarr[3, 1].xaxis.set_tick_params(labelsize=11)
    axarr[3, 1].yaxis.set_tick_params(labelsize=11)
    axarr[3, 1].set_ylim([-3,3.0])
    axarr[3, 1].yaxis.offsetText.set_visible(False)
    axarr[3, 1].text(-2.3, 1.8, 'Ti', fontsize=14, fontweight='bold')
    
    axarr[3, 2].scatter(yNi1, Ni1,c='r',s=100,marker='*')
    axarr[3, 2].scatter(yNi2, Ni2,c='b',s=100,marker='*')
    axarr[3, 2].scatter(yNi3, Ni3,c='g',s=100,marker='*')
    axarr[3, 2].set_xlim([-2.5,0.0])
    axarr[3, 2].xaxis.set_tick_params(labelsize=11)
    axarr[3, 2].yaxis.set_tick_params(labelsize=11)
    axarr[3, 2].set_ylim([-3,3.0])
    axarr[3, 2].yaxis.offsetText.set_visible(False)
    axarr[3, 2].text(-2.3, 1.8, 'Ni', fontsize=14, fontweight='bold')
    
    axarr[3, 3].scatter(yZn1, Zn1,c='g',s=100,marker='*')
    axarr[3, 3].scatter(yZn2, Zn2,c='g',s=100,marker='*')
    axarr[3, 3].scatter(yZn3, Zn3,c='g',s=100,marker='*')
    axarr[3, 3].set_xlim([-2.5,0.0])
    axarr[3, 3].xaxis.set_tick_params(labelsize=11)
    axarr[3, 3].yaxis.set_tick_params(labelsize=11)
    axarr[3, 3].set_ylim([-3,3.0])
    axarr[3, 3].yaxis.offsetText.set_visible(False)
    axarr[3, 3].text(-2.3, 1.8, 'Zn', fontsize=14, fontweight='bold')
    
    axarr[4, 0].scatter(yY, Y,c='g',s=100,marker='*')
    axarr[4, 0].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
    axarr[4, 0].set_ylabel("[X]",fontsize=10, fontweight='bold')
    axarr[4, 0].set_xlim([-2.5,0.0])
    axarr[4, 0].xaxis.set_tick_params(labelsize=11)
    axarr[4, 0].yaxis.set_tick_params(labelsize=11)
    axarr[4, 0].set_ylim([-3,3.0])
    axarr[4, 0].yaxis.offsetText.set_visible(False)
    axarr[4, 0].text(-2.3, 1.8, 'Y', fontsize=14, fontweight='bold')
    
    axarr[4, 1].scatter(yZr, Zr,c='g',s=100,marker='*')
    axarr[4, 1].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
    axarr[4, 1].set_xlim([-2.5,0.0])
    axarr[4, 1].xaxis.set_tick_params(labelsize=11)
    axarr[4, 1].yaxis.set_tick_params(labelsize=11)
    axarr[4, 1].set_ylim([-3,3.0])
    axarr[4, 1].yaxis.offsetText.set_visible(False)
    axarr[4, 1].text(-2.3, 1.8, 'Zr', fontsize=14, fontweight='bold')
    
    axarr[4, 2].scatter(yBa, Ba,c='g',s=100,marker='*')
    axarr[4, 2].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
    axarr[4, 2].set_xlim([-2.5,0.0])
    axarr[4, 2].xaxis.set_tick_params(labelsize=11)
    axarr[4, 2].yaxis.set_tick_params(labelsize=11)
    axarr[4, 2].set_ylim([-3,3.0])
    axarr[4, 2].yaxis.offsetText.set_visible(False)
    axarr[4, 2].text(-2.3, 1.8, 'Ba', fontsize=14, fontweight='bold')
    
    axarr[4, 3].scatter(yLa, La,c='g',s=100,marker='*')
    axarr[4, 3].set_xlabel("[Fe]",fontsize=10, fontweight='bold')
    axarr[4, 3].set_xlim([-2.5,0.0])
    axarr[4, 3].xaxis.set_tick_params(labelsize=11)
    axarr[4, 3].yaxis.set_tick_params(labelsize=11)
    axarr[4, 3].set_ylim([-3,3.0])
    axarr[4, 3].yaxis.offsetText.set_visible(False)
    axarr[4, 3].text(-2.3, 1.8, 'La', fontsize=14, fontweight='bold')
    
    for index, filename in enumerate(filenames):
        abunds,mus = historyabund(filename)
        h = profile2dict(filename)
        
        ind = np.argmin(h['log_Teff'])        
        
        keys = list(abunds.keys())
        tots = list(mus.keys())
        
        scaled = np.zeros(len(tots))
        
        for i,tot in enumerate(tots):
            if tot in keys:
                scaled[i] = np.log10(abunds[tot][ind])-np.log10(mus[tot])+12.15
                
        print('----------------------------------------------------------------')
        print('PARAMETERS OF INTEREST')
        print('----------------------------------------------------------------')
        print('C12/C13 ', abunds['c12'][ind]/abunds['c13'][ind])
        print('O16/O18 ', round(abunds['o16'][ind]/abunds['o18'][ind],2))
        print('C/O     ', round(abunds['C'][ind]/abunds['O'][ind],2))
        print('RCB Age ', round(h['star_age'][ind],2), 'years')
        print('Model   ', int(h['model_number'][ind]))
        
        tots_name = np.array(['Li','C','N','O',\
                              'F','Ne','Na','Mg',\
                              'Al','Si','S','Ca',\
                              'Sc','Ti','Ni','Zn',\
                              'Y','Zr','Ba','La',\
                                  'Fe'])
            
        tots_val = np.zeros(len(tots_name))-99
        tots_val[12] = 7.54
        
        for j,tot in enumerate(tots):
            for i,name in enumerate(tots_name):
                if name == tot and tot in keys:
                    tots_val[i] = scaled[j]
        
        Lodders = np.array([3.35,8.46,7.90,8.76,\
                            4.53,7.95,6.37,7.62,\
                            6.54,7.61,7.26,6.41,\
                            3.15,5.00,6.29,4.70,\
                            2.28,2.67,2.25,1.25,\
                                7.54])
        ele = tots_val - Lodders
        
        print('----------------------------------------------------------------')	
        print('SURFACE ABUNDANCES FOR SELECTED SPECIES (Solar Unit Abundances)')
        print('----------------------------------------------------------------')
        for idx,name in enumerate(tots_name):
            print(name,round(ele[idx],2))
            
        if not labels:
            labels.append(filename)
        
        axarr[0, 0].scatter(ele[20],ele[0],label=labels[index],s=100,marker='s')
        axarr[0, 1].scatter(ele[20],ele[1],s=100,marker='s')
        axarr[0, 2].scatter(ele[20],ele[2],s=100,marker='s')
        axarr[0, 3].scatter(ele[20],ele[3],s=100,marker='s')
        axarr[1, 0].scatter(ele[20],ele[4],s=100,marker='s')
        axarr[1, 1].scatter(ele[20],ele[5],s=100,marker='s')
        axarr[1, 2].scatter(ele[20],ele[6],s=100,marker='s')
        axarr[1, 3].scatter(ele[20],ele[7],s=100,marker='s')
        axarr[2, 0].scatter(ele[20],ele[8],s=100,marker='s')
        axarr[2, 1].scatter(ele[20],ele[9],s=100,marker='s')
        axarr[2, 2].scatter(ele[20],ele[10],s=100,marker='s')
        axarr[2, 3].scatter(ele[20],ele[11],s=100,marker='s')
        axarr[3, 0].scatter(ele[20],ele[12],s=100,marker='s')
        axarr[3, 1].scatter(ele[20],ele[13],s=100,marker='s')
        axarr[3, 2].scatter(ele[20],ele[14],s=100,marker='s')
        axarr[3, 3].scatter(ele[20],ele[15],s=100,marker='s')
        axarr[4, 0].scatter(ele[20],ele[16],s=100,marker='s')
        axarr[4, 1].scatter(ele[20],ele[17],s=100,marker='s')
        axarr[4, 2].scatter(ele[20],ele[18],s=100,marker='s')
        axarr[4, 3].scatter(ele[20],ele[19],s=100,marker='s')

    f.legend()
    f.align_ylabels(axarr[:,0])
    plt.show()
    
    
        
        
        
        
        

