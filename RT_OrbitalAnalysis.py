# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: July 2017
Last Update: 10 August 2017

@author: Olivia Hull

HOW THIS LIBRARY IS SET UP:
    FUNCTIONS:
    (1) get_pop(filename) will output:
            -a .csv file of the orbital populations with respect to time
            -a .txt file of the fourier transform of orbital populations vs. time
            -plots of the .csv and .txt files
        get_pop only requires the real-time TDDFT .log file as input, e.g. in the python terminal, call get_pop("myfilename.log") to run.
    (2) get_info(filename) will output:
        the filename, multiplicity, timestep, number of steps, and number of lines of orbitals in the file.
    *nothing in here will compute the S(w), dipole strength function. call RT_spect in beocat to do that.
    (3) get_orbnumber(filename) will output the number of lines of orbitals that are in the output file.
        for example, under "Time =", there is a list of all of the orbitals.
        get_orbnumber just gets the number of lines of orbitals that there are, which is necessary for other parts of the script.
    (4) get_pop_3(filename, spin, timestep, nsteps, orbcount) outputs a .csv file of the time and orbital occupation values,
        which is set up as [time, orbital1, orbital2, orbital3, ..., orbital n], with the number of rows corresponding to the number
        of time steps. get_pop_3 is ONLY called when the system has a spin multiplicity of 3.
    (5) get_pop_1(filename, spin, timestep, nsteps, orbcount) does the exact same thing as get_pop_3, except
        for systems which have a spin multiplicity of 1.
    (6) occufft(real_time_file) takes the .csv file of (time, orbital occupation numbers) and fourier transforms the data
        then outputs a .txt file of [frequency, intensity1, intensity2, ... , intensity n]
    (7) plot_data(file, HOMO, LUMO) will take a .csv file of the orbital populations OR a .txt file of the FT of orbital populations
        and then plot the data. The plotting parameters inside can be tweaked to plot more or less HOMOs/LUMOs
    (8) plot_triplet(filename, aHOMO, aLUMO, bHOMO, bLUMO) automates the plotting of the four triplet data files.

TO RUN:
    call:
        get_pop("yourfile.log")
        # NOTE: be sure to type "from orbital_pops_and_FT_v2 import *" into your python terminal first
    -First, get_pop calls get_info.
        -then a message will ask you if the displayed jbo information is correct. It will display the multiplicity, time step,
        number of steps, simulation length, and number of orbital lines in the file. If all these are correct, type 1 and press enter.
    -Then, get_pop will call get_pop_1 or get_pop_3 depending on the multiplicity value that it found.
    -Then, it will call plot_data for each of the .csv and .txt files generated.
    -If your system is a triplet, it will direct to plot_triplet, which will then call plot_data four times for the four
        separate triplet files output by get_pop_3 and occufft.

SUMMARY:
    get_pop('filename') to get orbital populations and FT of orbital populations and their plots
    dipole_overlay('filename') to get the dipole strength function (need x, y, z file to be in same directory)
"""

from __future__ import division
import csv
import numpy as np
from scipy.fftpack import fft, fftfreq
import matplotlib.pyplot as plt
#import matplotlib.axes as ax
import sys


def get_pop(filename):
    
    filename, spin, timestep, nsteps, orbcount = get_info(filename)
    
    if spin == 1:
        file_out, HOMO, LUMO = get_pop_1(filename,spin,timestep,nsteps, orbcount) # outputs csv file of RT orbital populations
        file_out_fft = occufft(file_out)
        plot_data(file_out, HOMO, LUMO)
        plot_data(file_out_fft, HOMO, LUMO)# outputs FT of RT orbital populations
        return HOMO, LUMO
    if spin == 3 or spin == 2:
        file_out_alpha, file_out_beta, aHOMO, aLUMO, bHOMO, bLUMO = get_pop_3(filename, spin, timestep, nsteps, orbcount)
        fout_alpha_fft = occufft(file_out_alpha)
        fout_beta_fft = occufft(file_out_beta)
        plot_triplet(filename, aHOMO, aLUMO, bHOMO, bLUMO)
#        fout_total_fft = occufft(file_out_total)
        return aHOMO, aLUMO, bHOMO, bLUMO
        
def get_info(filename):
    count = 0 # the third instance of "#" is the line we want
    spin_line = 0 # initializing? I think this line is superfluous
    info_line = [] # initialize the information line (third instance of #)
    with open(filename, 'r') as f:
        for line in f:
            if "Multiplicity =" in line:
                spin_line = line
            if "#" in line:
                count = count + 1
                if count == 3: # changed to 5 for new file format
                    while line[0:2] != ' -':
                        info_line.append(line)
                        line = next(f)
                    break
    
    spin_line = spin_line.split() 
    spin = spin_line[-1] # the spin value is the second to last of the line.
    spin=int(spin)
    
    # Below is just string processing
    for i in info_line:
        if "iop" in i:
            iop_find=i.split()
        
    iop_line = []
    for i in iop_find:
        if "5/" in i:
            iop_line.append(i)
    
    iop_info = []
    for i in iop_line:
        iop_info.append(i)
    
    iop_combine=' '
    for i in iop_info:
        iop_combine = iop_combine + i
    
    iop_info = iop_combine[5:-1]
    iop = iop_info.split(',')
    
    # the iop infos. See the iop printed sheet for what they mean
    for i in iop:
        if i[2:5] == '177':
            nsteps = float(i[6:])
        if i[2:5] == '134':
            time_as = float(i[7:])
        if i[2:5] == '144' and i[-2:] != '-1':
            orbs = float(i[6:])
            orbcount = orbs*2 + 2
            orbcount = int(np.ceil(orbcount / 6)) # six entries per row in log file
        if i[2:5] == '144' and i[-2:] == '-1':
            orbcount = get_orbnumber(filename)
            
    timestep = time_as/(1E5)
    sim_length = timestep*nsteps
    
    print('multiplicity= ', spin, 'timestep= ', timestep, 'fs ', 'no. of steps =', nsteps)
    print('This means the simulation length = ', int(np.ceil((sim_length))), 'fs')
    print('number of orbital lines = ', orbcount)
    print('')
    print('Please check that this is correct')
    
    #check = input('Enter 1 if yes, 0 if no: ')
    
    #if check == str(1):
    return filename, spin, timestep, nsteps, orbcount # pass these values into get_info, which directs them to get_pop_3 or get_pop_1, depending
    #if check == str(0):
    #   print('Error in function Get_Info. Call def get_pop_1(filename,spin, timestep, nsteps, count')
            

def get_orbnumber(filename):

    length = []
    counter = 0
    end_counter = 0

    with open(filename, 'r') as f:
        for line in f:
            if ("rbital occupation numbers:" in line):
                orbline = next(f)
                while end_counter == 0:# this stays true until end of file
                    if orbline[:] == '\n' or orbline[:] == '' or "Beta" in orbline:
                        end_counter += 1
                        break
                    else:
                        counter += 1
                        length.append(len(orbline.split()))
                        orbline = next(f)
                end_counter += 1

    
    # the entries in variable "length" tell of how many entries are in each line, so the smallest value will be the last line
    orbcount = len(length)
    return orbcount
                    
            # FIX THIS WITH SPLIT FUNCTION ... want to split the orbline, test if the type is integer.)
    
#    lastline = # want minimum value of the length string
                
def get_pop_3(filename, spin, timestep, nsteps, count):

    
    # alpha and beta _datalines are an array with the first index corresponding to the orbital line, and the second index corresponding to the time, e.g.
    # alpha_datalines[12][50000] prints out the orbital populations of the orbitals printed on the 12th orbital line of the output
    # at the 50,000th time step
    alpha_datalines = []
    beta_datalines = []

    for i in range(count):
        alpha_datalines.append([])
        beta_datalines.append([])
    
    # Get Alpha orbital info    
    with open(filename, 'r') as f:
        for line in f:
                if ("Alpha orbital occupation numbers:" in line):
                    alpha_val = next(f)
                    
                    while True:
                        if alpha_val[0] == '\n' or "Beta" in alpha_val:
                            break
                        if alpha_val[1] == ' ':
                            for i in range(count):
                                alpha_datalines[i].append(alpha_val)
                                alpha_val = next(f)
                        else:
                            break
    # Get beta orbital info
    with open(filename, 'r') as f:
        for line in f:
                if ("Beta orbital occupation numbers:" in line):
                    alpha_val = next(f)
                    
                    while True:
                        if alpha_val[0] == '\n':
                            break
                        try:
                            float_test = float(alpha_val.split()[0])
                        except:
                            e = sys.exc_info()[0]
                            break
                        if alpha_val[1] == ' ':
                            for i in range(count):
                                beta_datalines[i].append(alpha_val)
                                alpha_val = next(f)
                        else:
                            break
                        

    # Process/rearrange the data:
    steps = int(nsteps)
    
    time_list = [None for x in range(steps)]
    time_list[0] = float(0)
    
    time_list[1] = time_list[0] + timestep
    for i in range(2, steps):
        time_list[i] = time_list[i-1] + timestep
    
    alpha_all = []                    
    for j in range(steps):
        alpha_all.append(str(time_list[j]))
        for i in range(count):
            alpha_all[j] = alpha_all[j] + alpha_datalines[i][j]
    
    beta_all = []                    
    for j in range(steps):
        beta_all.append(str(time_list[j]))
        for i in range(count):
            beta_all[j] = beta_all[j] + beta_datalines[i][j]
    
    alpha_occ = []
    for i in range(steps):
        alpha_occ.append([float(x) for x in alpha_all[i].split() ])
    
    beta_occ = []
    for i in range(steps):
        beta_occ.append([float(x) for x in beta_all[i].split() ])
    
#    tot_occ = []
#    for i in range(steps):
#        tot_occ.append([])
#        tot_occ[i].append([])
#        tot_occ[i][0] = time_list[i]
#        for j in range(1, len(alpha_occ[0])):
#            tot_occ[i].append([])
#            tot_occ[i][j] = alpha_occ[i][j] + beta_occ[i][j]
    
    # Save the data; write to .csv files
    file_split= filename.split('.')
    file_out_alpha = file_split[0] + '_alpha' + '.csv'
    file_out_beta = file_split[0] + '_beta' + '.csv'
#    file_out_total = file_split[0] + '_total' + '.csv'
    
    with open(file_out_alpha, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(alpha_occ)
        
    with open(file_out_beta, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(beta_occ)
        
#    with open(file_out_total, 'w') as f:
#        writer = csv.writer(f)
#        writer.writerows(tot_occ)
                
    # Get the indices of the alpha HOMO LUMO and the beta HOMO LUMO            
    for i in range(len(alpha_occ[0])):
        if int(alpha_occ[0][i]) == 1 and int(alpha_occ[0][i+1]) == 0:
            aHOMO = i
            aLUMO = i + 1
            print('aHOMO = ', i, 'aLUMO = ', i+1)
            
    for i in range(len(beta_occ[0])):
        if int(beta_occ[0][i]) == 1 and int(beta_occ[0][i+1]) == 0:
            bHOMO = i
            bLUMO = i + 1
            print('bHOMO = ', i, 'bLUMO = ', i+1)

                
    return file_out_alpha, file_out_beta, aHOMO, aLUMO, bHOMO, bLUMO #, fileout_total
    
def get_pop_1(filename,spin, timestep, nsteps, count):
    
    alpha_datalines = []

    for i in range(count):
        alpha_datalines.append([])
        
    with open(filename, 'r') as f:
        for line in f:
                if ("Orbital occupation numbers:" in line):
                    alpha_val = next(f)
                    
                    while True:
                        if alpha_val[0] == '\n':
                            break
                        if alpha_val[1] == ' ':
                            for i in range(count):
                                alpha_datalines[i].append(alpha_val)
                                alpha_val = next(f)
                        else:
                            break
                        
    
    steps = int(nsteps)

    time_list = [None for x in range(steps)] # preallocate
    time_list[0] = float(0)
    
    time_list[1] = time_list[0] + timestep
    for i in range(2, steps):
        time_list[i] = time_list[i-1] + timestep
    
    alpha_all = []                    
    for j in range(steps):
        alpha_all.append(str(time_list[j]))
        for i in range(count):
            alpha_all[j] = alpha_all[j] + alpha_datalines[i][j]
            
    
    alpha_occ = []
    for i in range(steps):
        alpha_occ.append([float(x) for x in alpha_all[i].split() ])
        
    file_split= filename.split('.')
    file_out_alpha = file_split[0] + '.csv'
    
    
    with open(file_out_alpha, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(alpha_occ)
    
    for i in range(len(alpha_occ[0])):
        if int(alpha_occ[0][i]) == 2 and int(alpha_occ[0][i+1]) == 0:
            HOMO = i
            LUMO = i+1
            print('HOMO = ', i, ' LUMO = ', i+1)
        
    return file_out_alpha, HOMO, LUMO
    
    
def occufft(real_time_file):    

    rt = np.genfromtxt(real_time_file, delimiter=',')
    #rt=pd.read_csv(real_time_file, sep=',',header=None)

    t = rt[:,0]
    
    columns = len(rt[1,:])
    
    T = (t*41.341) # convert from au to fs

    for i in range(1,columns):
        orb_pop      = rt[:,i]
        # do the fourier transform 
        fw = fft(orb_pop)
    
        # determine frequency range
        n = len(fw)
        timestep = T[1] - T[0]              # spacing between time samples; assumes constant time step
        w = fftfreq(n,d=timestep)*2.0*np.pi # frequency list
        
        S = abs(fw)
        w = (w*27.2114)    # give frequencies in eV

        w = w[1:50000]
        S = S[1:50000]
        
        if i == 1:
            z = np.transpose([w,S])
        else:
            z = np.column_stack((z, S))
        
    f_out = real_time_file.split('.')
    f_out = f_out[0]
    f_out = f_out + '_occufft_data.txt'
    
    np.savetxt(f_out, z)
    return f_out

def plot_data(file, HOMO, LUMO):
    
    orbrange=10
    
    if '_z' in file:
        direction = 'z'
    if '_z' in file:
        direction = 'z'
    elif '_y' in file:
        direction = 'y'
    elif '_x' in file:
        direction = 'x'
    
    
    if file[-3:] == 'csv':
        data = np.genfromtxt(file, delimiter = ',')
    
    if file[-3:] == 'txt':
        data = np.genfromtxt(file)
    
   # len_data = len(data[0])
    
    # LUMO PLOTTING, population:
    plt.figure()
    ax = plt.subplot(111)
    # if you want all orbitals, to get all LUMOs, say
        # "for i in range(LUMO, len(data[0,:]))
    for i in range(LUMO,LUMO+orbrange): # previously : for i in range(LUMO,len_data):
        plt.plot(data[:,0], data[:,i], label='LUMO +' + str(i-LUMO), linewidth=0.5)
        plt.ylim([np.amin(data[0:60000,LUMO:]) - .05*np.amin(data[0:60000,LUMO:]),np.amax(data[0:60000,LUMO:]) + .05*np.amax(data[0:60000,LUMO:])])
        
    # Shrink current axis by 20%
    # previously ax.blah

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    fout = file[:-4]
    
    if file[-3:] == 'csv':
        plt.title('LUMO populations, %s-direction' %direction)
        plt.ylabel('Population')
        plt.xlabel('Time (fs)')
        plt.xlim([0,60])
        plt.show()
        savefile = fout + '_LUMO_pops.png'
        plt.savefig(savefile, format='png', dpi=300)
    
    if file[-3:] == 'txt':
        plt.xlim([0,8])
        plt.ylim([0,np.amax(data[:,1:]) + .01])
        plt.title('Fourier Transform of LUMO populations, %s-direction' %direction)
        plt.ylabel('Intensity')
        plt.xlabel('Energy (eV)')
        plt.show()
        savefile = fout + '_LUMO_FT.png'
        plt.savefig(savefile, format='png', dpi=300)
        #       plt.figure
    plt.show()
    # HOMO PLOTTING, population:
        
    plt.figure()
    ax = plt.subplot(111)
    for i in range(HOMO-orbrange,HOMO+1): # previously : for i in range(LUMO,len_data):
        plt.plot(data[0:60000,0], data[0:60000,i], label='HOMO - ' + str(HOMO-i), linewidth=0.5)
        #plt.ylim([np.amin(data[:,HOMO-5:HOMO]) - .00001*np.amin(data[:,HOMO-5:HOMO]), np.amax(data[:,HOMO-5:HOMO])])
    
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    if file[-3:] == 'csv':
        plt.show()
        plt.title('HOMO populations, %s-direction' %direction)
        plt.ylabel('Population')
        plt.xlabel('Time (fs)')
        plt.xlim([0,60])
        plt.show()
        savefile = fout + "_HOMO_pops.png"
        plt.savefig(savefile, format='png', dpi=300)
    
    if file[-3:] == 'txt':
        plt.xlim([0,8])
        plt.ylim([0,np.amax(data[:,1:]) + .1])
        plt.title('Fourier Transform of HOMO populations, %s-direction' %direction)
        plt.ylabel('Intensity')
        plt.xlabel('Energy (eV)')
        plt.show()
        savefile = fout + "_HOMO_FT.png"
        plt.savefig(savefile, format='png', dpi=300)
    
    # LUMO PLOTTING, FT:
        

def plot_triplet(filename, aHOMO, aLUMO, bHOMO, bLUMO):
    # supply the log file name
    
    f_prefix = filename[:-4]
    for i in ['_alpha.csv', '_beta.csv', '_alpha_occufft_data.txt', '_beta_occufft_data.txt']:
        file_in = f_prefix + i
        if 'alpha' in i:
            plot_data(file_in, aHOMO, aLUMO)
        if 'beta' in i:
            plot_data(file_in, bHOMO, bLUMO)
    
def time_dipole(filename):
    
    time_dipole = []
    
    with open(filename, 'r') as f:
        for line in f:
                if ("Time =" in line):
                    while True:
                        time_val = line
                        time_list = time_val.split()
                        time_fs = time_list[-4]
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
    
def cqRealTime(real_time_file,dipole_direction,kick_strength,damp_const):    
    '''
        (C) Joshua Goings 2016
            revised by Olivia Hull August 2017
        
        CQ_RealTime.py: a post-processing script for computing the absorption spectrum of
         Real Time Time Dependent SCF jobs in Chronus Quantum
        Computes the energy range, w (in eV) and dipole strength function S(w) for
         a given real time TD-SCF run. 
        real_time_file   ... type:string ; the RealTime_Dipole.csv file from a ChronusQ run
        dipole_direction ... type:char   ; which dipole moment contribution is computed (e.g. 'x','y', or 'z')
        kick_strength    ... type:float  ; in a.u., what was the applied field strength (e.g. 0.0001 au)
        damp_const       ... type:float  ; in a.u. of time, gives FWHM of 2/damp_const
    '''
    # for damp_const, want in units of eV. damping 27.211396 (1 eV in a.u.) --> 0.2 eV (use value 272.211396 for damping to get 0.2 eV FWHM --> this gives approx. experimental spectrum FWHM/variation)

    # in RT.com file, have field=x+$VAL, where $VAL converts to kick_strength of e.g. x+10 converts to 0.001 a.u.

    # this script takes dipole moments in the specified field direction and plots dipole moment vs time


    # chronusq file is CSV, also skip the header (first row)
    #    rt = np.genfromtxt(real_time_file,skip_header=1,delimiter=',')
    # rt = np.genfromtxt(real_time_file,delimiter=',')
    rt = np.genfromtxt(real_time_file,delimiter=',') # oah added. original is line directly above this one
    # choose which dipole axis you want
    if dipole_direction.lower() == 'x':
        direction = 1
    elif dipole_direction.lower() == 'y':
        direction = 2
    elif dipole_direction.lower() == 'z':
        direction = 3
    elif dipole_direction.lower() == 't':
        direction = 4
    else:
        print("Not a valid direction for the dipole! Try: x,y,z, or t for total")
        sys.exit(0)
    
    t      = rt[:,0] # put time values into its own array
    #print t
    # note 'z' is just generic dipole direction, converted from debye to au
    z      = rt[:,direction]*0.393456 # put dipole values into its own array, then multiply by conversion factor
    #   z      = rt[:,direction]*0.393456


    # scale dipole signal  
    zdif = (np.amax(z) - np.amin(z))/2.0 # np.amax(z) returns max dipole value in the list
    z = z - zdif # scale the dipole signal such that the middle value is at 0
    
    # add damping to give Lorenztian lineshape with FWHM of (2/damp_const)
    damp = np.exp(-(t-t[0])/damp_const)
    z = z * damp

    # pad signal with zeros. also not necessary, but makes spectra prettier
    #zero = np.linspace(0,0,10000)
    #z = np.hstack((z,zero))
   
    # do the fourier transform 
    fw = fft(z)
    
    # determine frequency range
    n = len(fw)                         # number samples, including padding
    timestep = t[1] - t[0]              # spacing between time samples; assumes constant time step
    w = fftfreq(n,d=timestep)*2.0*np.pi # frequency list
   
    fw_re = np.real(fw)                 # the real FFT frequencies
    #fw_im = (np.imag(fw))               # the imaginary FFT frequencies
    #fw_abs = abs(fw)                    # absolute value of frequencies # OAH (re^2 + im^2)^(1/2)w
    

    # 'correct' equation for dipole strength function assuming you did SCF in static field
    
    # below if/else statement added by OAH. If the first dipole term, z[0], is positive, then use the first equation. If the first dipole term is negative, then use the second equation
    if rt[0][direction] == abs(rt[0][direction]):
        S = (2.0*w*w*fw_re)/(3.0*np.pi*137*kick_strength)
    else:
        S = -(2.0*w*w*fw_re)/(3.0*np.pi*137*kick_strength)

    # 'correct' equation for dipole strength function assuming you did a delta kick
    #S = -(4.0*w*np.pi*fw_im)/(3.0*137*kick_strength)
    
    w = (w*27.2114)    # give frequencies in eV
    
    f_out_prefix = real_time_file[:-4]
    f_out = f_out_prefix + 'dip_stren.txt'
    np.savetxt(f_out, np.transpose([w,S]))
    
    return w, S, f_out

#if __name__ == '__main__':

#    Filename   = 'time_dip.csv', linewidth=0.5

#    xFilename   = 'h2o_x_RealTime_Dipole.csv'
#    yFilename   = 'h2o_y_RealTime_Dipole.csv'
#    zFilename   = 'h2o_z_RealTime_Dipole.csv'
    
#    kick        = 0.001 # depends on system
  # damping     = 272.11396 #OAH # anywhere between 50-250 usually works well
#    damping     = 272.11396
#    w, Sxx      = cqRealTime(xFilename,'x',kick,damping)
#    w, Syy      = cqRealTime(yFilename,'y',kick,damping)
#    w, Szz      = cqRealTime(zFilename,'z',kick,damping)

#    w, S      = cqRealTime(Filename,'o',kick,damping) # in function, takes 'o' as 2nd column rather than specify the actual direction, since the .csv file made with RT_spec.bash only harvests the requested dipole moment (so only two columns in the file)
#    np.savetxt('spectrum_ev_dipstren_test.txt', np.transpose([w,S]))  
    
#    plt.plot(w,Sxx+Syy+Szz,label='S')
# Can change these filenames and the title below, as well as x, y value ranges
##    plt.plot(w,S,label='S')
#    plt.ylim(-20,20)  # y range # OAH made range negative
#    plt.xlim(0,15)     # X range # OAH was (0.8)
#    plt.legend()
#    plt.show()
#    plt.savefig('quickplot_test.pdf')

def dipole_overlay(filename):
    # NEED TO CALL TIME_DIPOLE FUNCTION FIRST!!
    dirlist = ['_x', '_y', '_z']
    for index, item in enumerate(dirlist):
        next1 = index - 1
        next2 = index - 2
        if item in filename:
            filename1 = filename.replace(item, dirlist[next1])
            filename2 = filename.replace(item, dirlist[next2])
         
    damp_const = 272.11396
    kick_strength = 0.001
    
    #plt.figure()
    for direction in [filename, filename1, filename2]:
        if '_x_' in direction:
            time_v_dipolex = time_dipole(direction)
            w_x, S_x, f_out_x = cqRealTime(time_v_dipolex, 'x', kick_strength, damp_const)
            #plt.plot(w_x, S_x, label = 'x-direction')
        if '_y_' in direction:
            time_v_dipoley = time_dipole(direction)
            w_y, S_y, f_out_y = cqRealTime(time_v_dipoley, 'y', kick_strength, damp_const)
            #plt.plot(w_y, S_y, label = 'y-direction')
        if '_z_' in direction:
            time_v_dipolez = time_dipole(direction)
            w_z, S_z, f_out_z = cqRealTime(time_v_dipolez, 'z', kick_strength, damp_const)
            #plt.plot(w_z, -S_z, label = 'z-direction')
    plt.figure()
    plt.plot(w_x, S_x, label = 'x-direction')
    plt.plot(w_y, S_y, label = 'y-direction')
    plt.plot(w_z, S_z, label = 'z-direction')
    
    w = w_x # redefinig for clarity. w should not change between w_x, w_y, w_z regardless
    S_total = S_x + S_y + S_z
    plt.plot(w, S_total, label = 'Total')
    
    plt.ylim([0,np.amax(S_total[0:500]) + .5])
    plt.xlim([0,7])
    plt.legend(loc='upper right')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Dipole Strength')
    plt.title('Dipole Strength Function')
    plt.show()
    savefile = filename[:-4] + "_dipole_overlay.eps"
    plt.savefig(savefile, format='eps', dpi=1000)
    
def max_peaks(file):
    
    if file[-3:] == 'csv':
        data = np.genfromtxt(file, delimiter = ',')
    
    if file[-3:] == 'txt':
        data = np.genfromtxt(file)
    
    maxs = np.amax(file[:])
    [i for i, j in enumerate(a) if j == maxs]
    

def occufft_actual_time(real_time_file):    

    rt = np.genfromtxt(real_time_file, delimiter=',')
    #rt=pd.read_csv(real_time_file, sep=',',header=None)

    t = rt[:,0]
    
    columns = len(rt[1,:])
    
    T = (t*41.341) # convert from au to fs

    for i in range(1,columns):
        orb_pop      = rt[:,i]
        # do the fourier transform 
        fw = fft(orb_pop)
    
        # determine frequency range
        n = len(fw)
        timestep = T[1] - T[0]              # spacing between time samples; assumes constant time step
        w = fftfreq(n,d=timestep)*2.0*np.pi # frequency list
        
        S = abs(fw)
        w = (w*27.2114)    # give frequencies in eV

        w = w[1:50000]
        S = S[1:50000]
        
        if i == 1:
            z = np.transpose([w,S])
        else:
            z = np.column_stack((z, S))
        
    f_out = real_time_file.split('.')
    f_out = f_out[0]
    f_out = f_out + '_occufft_data.txt'
    
    np.savetxt(f_out, z)
    return f_out

def Autorunner(filename):
    a = 0
    return a    