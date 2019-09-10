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
import inspect
import importlib
import os

'''
Things to do:
    6. Add Ag-N distance function (just track one N atom rather than COM or anything like that)
'''


def my_main(filename, eF_type):
    pos1=9
    pos2=10
    a1='N'
    a2='N'
    fname = filename[:-4] # get rid of .log
    coord_lines = get_nuc_coords(filename)
    td = getDistance(coord_lines, pos1, pos2)
    params = fname.split("_")
   # eF_type = params[2][0:4] # e500, etc.
    orient = params[0] # 'par,' etc.
    polar = params[1] # e-field polarization, long or trans

    eFields= {
      "e500" : "0.05",
      "e450" : "0.045",
      "e400" : "0.04",
      "e300" : "0.03",
      "e200" : "0.02",
      "e100" : "0.01",
      "e10"  : "0.001"
      }

    orientations = {
      "par" : "parallel",
      "perp" : "perpendicular"
      }

    polarizations = {
      "long" : "longitudinal",
      "trans" : "transverse"
      }
    
    
    efield_label = eFields[eF_type]
    orientation_label = orientations[orient]
    polarization_label = polarizations[polar]
    calc_params = efield_label + " a.u. " + polarization_label+"-polarized electric field, " \
                  + orientation_label + " orientation"

    figOut=fname + "_" + eF_type + ".png"
    plotDistance(td, a1, a2, figOut, calc_params)


def get_nuc_coords(filename):
    ''' Pull the nuclear coordinates from an ehrenfest dynamics calculation in Gaussian
    input: filename (string)
    output: coordinate file, csv format.
    Name of coordinate file is "filename_nc.csv" '''
    
    count = 1
    coord_lines = [count]
    
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
    f.close()
    
    return coord_lines                    
    
    '''Now all of the coordinate data is in coord_lines'''

    
def getDistance(coord_lines, pos1, pos2):
    #pos1 = 9 # 9th position in coord list
    #pos2 = 10 # 10th position in coord list
    
    # find number of atoms
    nA = 0
    while type(coord_lines[nA+1]) == list:
        nA = nA+1
    
    # preallocate
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
        t_fs = t_fs + 0.0005999 # put time in units of fs  

    td = np.zeros((len(t), 2))
    td[:, 0] = np.array(t)
    td[:, 1] = np.array(d)

    return td

def plotDistance(td, a1, a2, figOut, calc_params):
    '''
    '''
    a1 = 'N'
    a2 = 'N'
    
    plt.figure()
    
    plt.plot(td[:,0], td[:,1])
    plt.xlabel("Time (fs)")
    plt.ylabel("Distance (Angstrom)")

    plt.suptitle(a1+"-"+a2+" bond distance over time", y=1.05, fontsize=18)
    plt.title(calc_params, fontsize=10)    

    plt.savefig(figOut, format='png', dpi=300, bbox_inches='tight')
    

filename = sys.argv[1]
eF_type = sys.argv[2]

my_main(filename, eF_type)
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
    
    

    
