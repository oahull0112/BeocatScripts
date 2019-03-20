#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: July 2018
@author: Olivia Hull

HOW THIS LIBRARY IS SET UP:
    FUNCTIONS:
        (1) RT_offdiag_control(filename, coupling_array)
                Controls which functions are called on the data.
                Use this function to run analysis.
                        INPUT:
                        -filename: input the filename as a string, e.g. 'Ag6.log'
                        -coupling_array: input the couplings of interest as an array
                        e.g. RT_offdiag_control('Ag6.log', [[64, 65], [64,66]])
        (2) time_dipole(filename):
                Grabs the time and dipole data from the file, then outputs as a .csv file
        (3) plot_time_dipole(dipole_data):
                Plots the dipole data and saves as a png
        (4) RT_offdiag(filename):
                Grabs ALL off-diagonal data and inputs into a python dump file
        (5) plot_RT_offdiag:
                Plots the off-diagonals of interest, which were specified in the variable coupling_array
        (6) FT_offdiags:
                Computes the fourier transform of the off-diagonals specified in coupling_array
        (7) plot_FT_offdiags:
                Plots the fourier transform of the couplings specified in coupling_array
    TO RUN:
        Call an interactive job session via the terminal command:
                srun -p killable.q --constraint=avx2 --nodes=1 ==ntasks-per-node=1 --mem-per-cpu=4G ==pty /bin/bash
                NOTE:   You can specify more memory if necessary, but the code is not parallelized so calling more than 1 processor will do nothing.
        Then type "python2"
        Inside python2, type "from OffDiag import *"
        then input the command "RT_offdiag_control(filename, coupling_array), where you have specified the filename as a string and the coupling array as an array.
                e.g. RT_offdiag_control('Ag6.log', [[64, 65], [64,66]]) 
   PLOTTING:
        If you do not like the look of a plot, e.g. you need to change the x or y axis values, legend location, etc.
        Then you must do that inside the plotting function itself by going into OffDiag.py and changing the value.
"""


from __future__ import division
import csv
import numpy as np
from scipy.fftpack import fft, fftfreq
import matplotlib.pyplot as plt
#import matplotlib.axes as ax
import sys
import json



def RT_offdiag_control(filename, coupling_array):
    ''' "control center" for the off-diagonal type calculation analysis'''
    
    fout_RT_plot = filename[:-4] + '_RT_Offdiag.png' # you can change this to whatever; the output plot filename for couplings plot
    fout_FT_plot = filename[:-4] + '_FT_Offdiag.png'
    
    dipole_data = time_dipole(filename) # outputs the file name containing the dipole information
    plot_time_dipole(dipole_data) # makes the dipole plot
    coupling_data, coupling_data_file = RT_offdiag(filename) # outputs the actual information matrix and the file name containing coupling info
    data_array, time_data, find_indices = plot_RT_offdiag(coupling_data, dipole_data, coupling_array, fout_RT_plot)
    FT_out, FT_data = FT_offdiags(filename, data_array, time_data, find_indices)
    plot_FT_offdiags(FT_data, coupling_array)
    
def time_dipole(filename):
    '''working'''
    time_dipole = []
    
    with open(filename, 'r') as f:
        for line in f:
                if ("Time =" in line):
                    while True:
                        time_val = line
                        time_list = time_val.split()
                        time_fs = time_list[-2]
                        next(f)
                        dipole_line = next(f)
                        dip_list = dipole_line.split()
                        dip_vals = dip_list[1::2] # CHECK THIS!
                        dip_vals.insert(0, time_fs)
                        time_dipole.append(dip_vals) # list of dipole values in order [X, Y, Z, TOTAL]
                        break
    
    time_dipole.pop(0)
    #    del time_dipole[-1]
    file_out = filename[:-4] + '_time_dipole.csv'
    
    with open(file_out, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(time_dipole)

    return file_out





def plot_time_dipole(dipole_data):
    '''working'''
    plot_data = np.genfromtxt(dipole_data, delimiter = ',')
    j=1
    for i in ['x', 'y', 'z']:
        plt.plot(plot_data[:,0], plot_data[:,j], label=i)
        j=j+1
    plt.legend()
    plt.xlabel('time (a.u.)')
    plt.ylabel('dipole')
    plt.title('Dipole Moment')
    
    


def RT_offdiag(filename):
    '''
    working; Pull out the off-diagonal element information from the output file, return as a .csv
    data is set up as data[index1][index2] where index1 is the time step and index 2 is the off-diagonal pair
    '''
    
    data = []
    counter = -1

    with open(filename, 'r') as f:
        for line in f:
            if ("MO Density OV block:" in line):
                counter = counter + 1
                coupling_line = next(f)
                data.append([])
                while True:
                    if coupling_line[0:21] == '                     ':
                        split_line = coupling_line.split()
                        #data[counter].append(coupling_line.split())
                        data[counter].append([float(x) for x in coupling_line.split()])
                        coupling_line = next(f)
                    else:
                        break
            
    
    file_out = filename[:-4] + '_off_diags.list'
    
    with open(file_out, 'w') as out_file:
        json.dump(data, out_file) # to pull this back up, use:
        #with open(file_out, 'r') as in_file:
            #your_list = json.load(in_file)
                
    return data, file_out





def RT_offdiag_limited(filename, coupling_array):
    '''just pull out the couplings of interest'''
    data = []
    counter = -1

    with open(filename, 'r') as f:
        for line in f:
            if ("MO Density OV block:" in line):
                counter = counter + 1
                coupling_line = next(f)
                data.append([])
                while True:
                    if coupling_line[0:21] == '                     ':
                        split_line = coupling_line.split()
                        #data[counter].append(coupling_line.split())
                        data[counter].append([float(x) for x in coupling_line.split()])
                        coupling_line = next(f)
                    else:
                        break
            
    
    file_out = filename[:-4] + '_off_diags.list'
    
    with open(file_out, 'w') as out_file:
        json.dump(data, out_file) # to pull this back up, use:
        #with open(file_out, 'r') as in_file:
            #your_list = json.load(in_file)
                
    return data, file_out
    




def plot_RT_offdiag(data, dipole_data, coupling_array, fout):
    ''' Plotting for the RT off diag values. 
    data: input as either the raw python data list or as the out put filename of the data.
        It will be faster to input as the raw python data, if it is already in memory from the RT_offdiag function
        Otherwise, use the output file name.
    coupling_array: Provide the couplings you want to see in an array. For example, 
        coupling_array = [[1,2],[1,3],[1,4]] would plot the couplings between 
        (1,2), (1,3), (1,4) on the same plot
    fout: give the name of the saved output file you want
    '''
    
    if type(data) == list:
        True
    elif type(data) == str:
        with open(data, 'r') as in_file:
            data = json.load(in_file)
    
    data_array = np.array(data)
    
    find_indices = []
    #coupling_array = [[64, 65], [64, 66], [64, 67], [64, 68], [64, 69], [64, 70], [64, 71], [64, 72], [64, 73], [64, 74], [64, 75],[64,76]]
    
    time_data = np.genfromtxt(dipole_data, delimiter = ',') # dipole data is dipole file name
    time_data = time_data[:,0]
    
    # this finds the indices of the couplings of interest
    # Since the couplings may not start at orbital 1, have to correct for this in the array
    for pair in coupling_array:
        counter = 0
        for coupling in data[0]:
            if int(coupling[0]) == int(pair[0]) and int(coupling[1]) == int(pair[1]):
                find_indices.append(counter)
            counter = counter + 1
    
    plot_counter=0
    plt.figure()
    for pair in coupling_array:
        plt.plot(time_data[0:-1], data_array[:, find_indices[plot_counter], 2], label=str(pair[0]) + ', ' +str(pair[1]))
        plot_counter = plot_counter + 1

    
    plt.title('Off-diagonal couplings')
    plt.ylabel('Coupling')
    plt.xlabel('Step')
    plt.legend()
    plt.savefig(fout, format='png', dpi=300, bbox_inches='tight')
    
    return data_array, time_data, find_indices
    




    
def FT_offdiags(filename, data_array, time_data, find_indices):
    '''Calculate the fourier transform of the off-diagonals'''

    t = time_data    
    T = (t*41.341) # convert fs to au

    plot_counter = 0
    for pair in find_indices:
        coupling = data_array[:, find_indices[plot_counter], 2]
        # do the fourier transform 
        fw = fft(coupling)
    
        # determine frequency range
        n = len(fw)
        timestep = T[1] - T[0]              # spacing between time samples; assumes constant time step
        w = fftfreq(n,d=timestep)*2.0*np.pi # frequency list
        
        S = abs(fw)
        w = (w*27.2114)    # give frequencies in eV

        w = w[1:50000]
        S = S[1:50000]
        
        if plot_counter == 0:
            FT_data = np.transpose([w,S])
        else:
            FT_data = np.column_stack((FT_data, S))
            
        plot_counter = plot_counter+1
        
    f_out = filename.split('.')
    f_out = f_out[0]
    f_out = f_out + '_fft_data.txt'
    
    np.savetxt(f_out, FT_data)
    return f_out, FT_data





def plot_FT_offdiags(FT_data, coupling_array):
    
    if type(FT_data) == str:
        FT_data = np.genfromtxt(FT_data)

    plt.figure()
    ax = plt.subplot(111)

    for pair in range(len(coupling_array)): # previously : for i in range(LUMO,len_data):
        plt.plot(FT_data[:,0], FT_data[:,pair+1], label=str(coupling_array[pair][0]) + ", " + str(coupling_array[pair][1]), linewidth=0.5)
        
    #plt.ylim([np.amin(data[0:60000,LUMO:]) - .05*np.amin(data[0:60000,LUMO:]),np.amax(data[0:60000,LUMO:]) + .05*np.amax(data[0:60000,LUMO:])])

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylim([0, 10])
    plt.xlim([0, 10])
    
    plt.savefig('FT_plot.png', format='png', dpi=300, bbox_inches='tight')
    #fout = file[:-4]
    




def generic_plot(data, columns, labels, ylimit, xlimit, title):
    '''input columns as an array, i.e. [0, 2,5,6] would plot col. 0 against 2, then col 0 against 5, then col 0 against 6
    labels=["65, 66", "66, 67"]
    title="THis is my title"
    '''
    
    if type(data) == str:
        data = np.genfromtxt(data)
        
    plt.figure()
    #ax = plt.subplot(111)
    
    counter=0
    for yval in columns[1:]:
        plt.plot(data[:, columns[0]], data[:, yval], label=str(labels[counter]))
        counter = counter + 1
    
    plt.legend()
    plt.xlabel('Frequency')
    plt.ylabel('Intensity')
    plt.title(title)
    plt.ylim(ylimit)
    plt.xlim(xlimit)





def grab_time_data(filename):
    return 0
    



def combine_restarted_files(filename1, filename2, outputfile):  
# Get file inputs
# Example:  firstfile.log secondfile.log combinedoutput.log

#first, second, output = sys.argv[1], sys.argv[2], sys.argv[3]
    first = filename1
    second = filename2
    output = outputfile
    # Set expression to search
    #exp = 'Electronic Dynamics'
    exp = 'Step:     66245'
    # Get last expression of file 1
    f = open(first)
    for line_number, line in enumerate(f):
        if(exp in line):
            last_step = (line_number,line)
    
    # Find matching expression in file 2
    s = open(second)
    for line_number, line in enumerate(s):
        if(last_step[1] in line):
            match = (line_number,line)
            break
    
    # Obtain correct pieces of each file and combine lines
    f = open(first)
    lines = f.readlines()
    lines = lines[:last_step[0]]
    
    s = open(second)
    s_lines = s.readlines()
    s_lines = s_lines[match[0]:]
    lines.extend(s_lines)
    
    # Open output file to write
    o = open(output,'w')
    for item in lines:
        o.write('%s' % item)
        
def combine_restarts(first, second, output):
    '''
    '''
