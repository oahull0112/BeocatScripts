# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 17:09:37 2019

@author: olivia
"""


from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from scipy.fftpack import fft, fftfreq
from scipy.signal import find_peaks
#import inspect
#import importlib
#import os

'''
Things to do:
    6. Add Ag-N distance function (just track one N atom rather than COM or anything like that)
'''


#def my_main(filename, eF_type, pos1, pos2, a1, a2, doFT):
def my_main(filename, pos1, pos2, a1, a2, doFT):
    fname = filename[:-4] # get rid of .log
    coord_lines, E_cons, eInfo_list = get_nuc_coords(filename)
    eInfo = parseEInfoList(eInfo_list)
    eInfo = str(eInfo)
    print("eInfo: ", eInfo)
    td, avgBL = getDistance(coord_lines, pos1, pos2)
    calc_params, figOut = getCalcParams(fname, eInfo)
    plotDistance(td, a1, a2, figOut, calc_params, avgBL)
    if doFT == True:
        FT_data = getDistFT(td)
        plotDistFT(FT_data, a1, a2, figOut, calc_params)
        plotECons(E_cons, figOut, calc_params) # only need to plot ECons once

def getCalcParams(fname, eInfo):
#pulls calculation parameters for plotting purposes
    params = fname.split("_")
   # eF_type = params[2][0:4] # e500, etc.
    orient = params[0] # 'par,' etc.
    polar = params[1] # e-field polarization, long or trans

    eFields= {
       "0.05" : "e500", 
       "0.045": "e450",
       "0.04" : "e400",
       "0.03" : "e300",
       "0.02" : "e200",
       "0.01" : "e100",
       "0.001": "e10" 
      }

    orientations = {
      "par" : "parallel",
      "perp" : "perpendicular"
      }

    polarizations = {
      "long" : "longitudinal",
      "trans" : "transverse"
      }
    
    eF_type = eFields[eInfo]
    efield_label = eInfo
    orientation_label = orientations[orient]
    polarization_label = polarizations[polar]
    calc_params = efield_label + " a.u. " + polarization_label+"-polarized electric field, " \
                  + orientation_label + " orientation"

    figOut=fname + "_" + eF_type + ".png"

    return calc_params, figOut

def get_nuc_coords(filename):
    ''' Pull the nuclear coordinates from an ehrenfest dynamics calculation in Gaussian
    input: filename (string)
    output: coordinate file, csv format.
    Name of coordinate file is "filename_nc.csv" '''
    
    count = 1
    coord_lines = [count]
    E_cons = []
    eInfo_list = []
    
    with open(filename, 'r') as f:
        for line in f:
            if ("Number     Number       Type             X           Y           Z" in line):
                next(f)
                coord = next(f)
                while True:
                    if coord[1] == '-': # signals end of coordinate lines
                        count += 1
                        coord_lines.append(count)
                        break
                    else:
                        coord_lines.append([coord.split()[1], coord.split()[3], coord.split()[4], coord.split()[5]])
                        coord = next(f)
            if (' Energy = ' in line ):
                E_cons.append(float(line.split()[-1]))
            if ('                            Ex =' in line):
                e_line = line
                for i in range(3):
                    eInfo_list.append(e_line.split()[2])
                    e_line = next(f)

    f.close()
    
    return coord_lines, E_cons, eInfo_list                    
    
    '''Now all of the coordinate data is in coord_lines'''

def parseEInfoList(eInfo_list):
    for eF in eInfo_list:
        if float(eF) != 0:
            eInfo = float(eF)
    return eInfo
    
def getDistance(coord_lines, pos1, pos2):
    #pos1 = 9 # 9th position in coord list
    #pos2 = 10 # 10th position in coord list
    
    # find number of atoms
    nA = 0
    while type(coord_lines[nA+1]) == list:
        nA = nA+1
    
    x1 = [] ; x2 = [] ; xd = []
    y1 = [] ; y2 = [] ; yd = []
    z1 = [] ; z2 = [] ; zd = []
    
    d = []
    t = []
    
    #maxiter = int(np.floor(len(coord_lines)/i))
    
    # get coordinates, compute distance at each time step
    ti = 0 # "time index"
    t_fs = 0
    j = 0
    i = nA + 1 # indexer
    while j <= len(coord_lines)-i:
        x1.append(float(coord_lines[pos1 + j][1]))
        y1.append(float(coord_lines[pos1 + j][2]))
        z1.append(float(coord_lines[pos1 + j][3]))
        
        x2.append(float(coord_lines[pos2 + j][1]))
        y2.append(float(coord_lines[pos2 + j][2]))
        z2.append(float(coord_lines[pos2 + j][3]))
        
        xd.append((x1[ti] - x2[ti])**2)
        yd.append((y1[ti] - y2[ti])**2)
        zd.append((z1[ti] - z2[ti])**2)
        
        d.append(np.sqrt(xd[ti] + yd[ti] + zd[ti]))
        t.append(t_fs)
        
        j =  j + i
        ti = ti + 1 # time step index  
        t_fs = t_fs + .005 # put time in units of fs  
        # NOTE: this is currently a hacky bandaid to get fs
        # if your nuclear time step (stepsize in input file) is not
        # stepsize=1000, then you must manually change this value

    td = np.zeros((len(t), 2))
    td[:, 0] = np.array(t)
    td[:, 1] = np.array(d)
    avgBL = np.mean(d[4010:]) # cut out the steps that have the field on

    return td, avgBL

def plotDistance(td, a1, a2, figOut, calc_params, avgBL):
    '''
    https://stackoverflow.com/questions/1713335/peak-finding-algorithm-for-python-scipy
    '''
    
    plt.figure()
    
    plt.plot(td[:,0], td[:,1])
    plt.xlabel("Time (fs)")
    plt.ylabel("Distance (Angstrom)")

    peaks, _ = find_peaks(td[:,1], prominence=.01)
    difs = []
    for i in range(len(peaks)-1):
        difs.append(td[:,0][peaks[i+1]]-td[:,0][peaks[i]])
    avgBL =  "%.3f" % avgBL
    plt.suptitle(a1+"-"+a2+" distance over time (avg: "+ avgBL +" Angstrom)" , y=1.05, fontsize=18)
    plt.title(calc_params, fontsize=10)    

    plt.savefig(figOut[:-4]+"_"+a1+a2+".png", format='png', dpi=300, bbox_inches='tight')

def getDistFT(td):
    '''
    td is already in units of fs
    '''
    td = td[4010:] # cut out the electric field part
    fw = fft(td[:,1])
    n = len(td[:,1])
    t = td[:,0]
    timestep = t[0] - t[1]
    w = fftfreq(n, d=timestep) #*2.0*np.pi
    w = w*(1e15)/(2.9979e10)
    S = abs(fw)
# convert this to cm-1:
# also, put in the average bond length
# also make plots for Ag-N
# get energy conservation as well
   # w = (w*27.2114) # give freqs in eV

    FT_data = np.transpose([w,S])
    #FT_data = [w, S]

    return FT_data

def plotDistFT(FT_data, a1, a2, figOut, calc_params):
    plt.figure()
    pk=[]
    peaks, _ = find_peaks(FT_data[:,1], prominence=10)
    yvals=FT_data[:,1]
    xvals=FT_data[:,0]
    for peak in peaks:
        if peak > 500:
            pk.append(peak)
    mpk=len(xvals) - 1
    for peak in pk:
        if yvals[peak] >= yvals[mpk] and xvals[peak] > 0:
            mpk = peak
    ax = plt.subplot(111)
    plt.plot(xvals, yvals)
    plt.xlim([0,5000])
    plt.ylim([0, yvals[mpk]+10])
    plt.xlabel("Wavenumber (cm-1)")
    plt.ylabel("Intensity (a.u.)")
    plt.suptitle("Fourier transform of "+ a1+"-"+a2+" bond distance over time", y=1.05, fontsize=18)
    label_peak = "%.3f" % xvals[mpk]
    plt.title(calc_params +", peak: "+ label_peak + " cm-1", fontsize=10)
    plt.savefig(figOut[:-4]+"_"+a1+a2+"_FT.png", format='png', dpi=300, bbox_inches='tight')

def plotECons(E_cons, figOut, calc_params):
    plt.figure()
    ax = plt.gca()
#    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    plt.plot(E_cons)
    plt.gca().yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:.3f}"))
    plt.suptitle("Energy conservation", y=1.05, fontsize=18)
    plt.ylabel("Energy (a.u.)")
    plt.xlabel("Timestep")
    plt.title(calc_params)
    plt.savefig(figOut[:-4]+"_ECons.png", format='png', dpi=300, bbox_inches='tight')



filename = sys.argv[1]
pos1=9;pos2=10;a1='N';a2='N';doFT=True
my_main(filename, pos1, pos2, a1, a2, doFT)
pos1=8;pos2=9;a1='Ag';a2='N'; doFT=False
my_main(filename, pos1, pos2, a1, a2, doFT)
#!/usr/bin/env python
#make executable in bash chmod +x PyRun
# https://stackoverflow.com/questions/3987041/run-function-from-the-command-line

#if __name__ == "__main__":
#    cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
#    if cmd_folder not in sys.path:
#        sys.path.insert(0, cmd_folder)
#
#    # get the second argument from the command line      
#    methodname = sys.argv[1]
#
#    # split this into module, class and function name
#    modulename, classname, funcname = methodname.split(".")
#
#    # get pointers to the objects based on the string names
#    themodule = importlib.import_module(modulename)
#    theclass = getattr(themodule, classname)
#    thefunc = getattr(theclass, funcname)
#
#    # pass all the parameters from the third until the end of 
#    # what the function needs & ignore the rest
#    args = inspect.getargspec(thefunc)
#    z = len(args[0]) + 2
#    params=sys.argv[2:z]
#    thefunc(*params)
    
    

    
