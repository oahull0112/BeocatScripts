#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 10:59:59 2018

@author: oah
"""


from __future__ import division
import csv
import numpy as np
from scipy.fftpack import fft, fftfreq
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
#import matplotlib.axes as ax
import sys
from mpl_toolkits import mplot3d


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
                        
    
    '''Now all of the coordinate data is in coord_lines'''
    
def view_coords(coord_lines):
    
    '''PUT LINES TO LOAD IN DATA HERE
    coord_lines = 
    '''
    ax = plt.axes(projection='3d')
    
    
    x_coords = []
    y_coords = []
    z_coords = []
    i = 0
    while i < len(coord_lines):
        if type(coord_lines[i]) == int:
            i = i+1
            '''PLOT STUFF'''
            #ax.scatter3D(xdata, ydata,zdata, c=zdata, cmap='Greens')
            #plt.cla()
            ax.set_xlim3d(-10, 10)
            ax.set_ylim3d(-10,10)
            ax.set_zlim3d(-10,10)
            ax.scatter3D(x_coords, y_coords, z_coords, depthshade=False)
            plt.title('step = ' + str(i))
            plt.pause(.001)
            x_coords = []; y_coords = []; z_coords = []
        else:
            x_coords.append(float(coord_lines[i][1]))
            y_coords.append(float(coord_lines[i][2]))
            z_coords.append(float(coord_lines[i][3]))
            i = i + 1
