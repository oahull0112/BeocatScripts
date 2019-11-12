#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 16:24:43 2019
@author: Olivia Hull

Calculate the percent composition of a molecular orbital of its atomic orbitals.
This package is written restrictively for Ag6-N2. Any deviation from my standard Ag6-N2 coordinate input
will need to be accounted for in the package. See relevant function comments

NOTE: In order to run this, your Gaussian calculation must be run with the pop=full keyword
      in order for it to actually print out all of the MO coefficient information
      
WARNING: this code is kind of garbage. Need to manually change in TWO functions if not using N9, N10 (for Ag8N2)
         originally this had to be changed from N7, N8 for Ag6N2 and I haven't changed it yet to make everything
         more general for any system. To fix the code, need to track down the bugs in the get_raw_coefs routine.
         It would also be possible to make a wrapper function that just runs one "chunk" at a time automatically
         and then concatenates the results.
WARNING: See the warning in find_percent_composition!!! Run one "chunk" of MOs at a time or it will fail.
"""

import numpy as np

def find_percent_composition(filename, MO_array):
    '''run this function to get the total result. Outputs the entire MO coefficient array 'data'
    and the partitions into N7 and N8.
    input:
        -filename: name of the Gaussian output file as a string
        -MO_array: array of the MO indices that you want the compositions of, e.g. MO_array = [1, 2, 3] would give
                   the composition of the first three orbitals.
    
    Note: this code is buggy. It will fail if you specify MO coefs in the MO_array if those MO coefs
    are not in the same "chunk" of the Gaussian output file
    To get around this, just make multiple calls...Run lines = determine_indexing_lines(MO_array)
    and if lines gives multiple inner arrays, run individual calculations for each inner array
    '''
    
    MO_array = sorted(MO_array)
    line_count = get_line_count(filename)
    raw_MO_data = get_raw_coefs(MO_array, line_count, filename)
    clean_array = make_initial_array(raw_MO_data, line_count)
    lines = determine_indexing_lines(MO_array)
    data_array = cleanup_raw_coefs(raw_MO_data, line_count, MO_array, clean_array, lines)
    N7, N8 = N2_CSPA_contributions(data_array, MO_array, line_count)
    
    #return data_array, N7, N8
    return N7, N8

def get_line_count(filename):
    ''' determines number of AOs (= number of MOs)'''
    
    line_count = 0
    
    with open(filename,'r') as f:
        for line in f:
            if ("     Eigenvalues -- " in line):
                MO_line = next(f)
                while True:
                    if MO_line[0:19] != '                   ':
                        line_count +=1
                        MO_line = next(f)
                    else:
                        break
                break
    
    return line_count

def get_raw_coefs(MO_array, line_count, filename):
    ''' pulls the raw coef data from gaussian output file and outputs the raw lines'''
    data = []
    MO = 0
    orb_line = []
    
    with open(filename, 'r') as f:
        for line in f:
            if ("Molecular Orbital Coefficients:" in line):
                coupling_line = next(f)
                while MO < len(MO_array):
                    if str(MO_array[MO]) in coupling_line.split() and str(MO_array[MO]) not in orb_line:
                        for i in coupling_line.split():
                            orb_line.append(i)
                        for i in range(line_count + 3):
                            data.append(coupling_line)
                            coupling_line = next(f)
                        MO += 1
                    if str(MO_array[MO]) in orb_line:
                        MO += 1
                        continue
                    else:
                        for i in range(line_count + 3):
                            coupling_line = next(f)

    return data

def locate_MO_index(MO, raw_MO_data_line):
    '''Quick function used to pull the correct data points for the cleanup function'''
    
    counter = 0
    split_MO_line = raw_MO_data_line.split()
    
    for orb in split_MO_line:
        if str(MO) == orb:
            MO_col = counter
        counter = counter+1
        
    return MO_col # indexes from 0 (MO_col = 2 means third column)

def make_initial_array(raw_MO_data, line_count):
    '''Sets up the cleaned-up array without adding the MO data yet. Just grabs the first four columns, which
    contain info on the AO function names and atoms'''
    
    clean_array = []
    
    for val in range(line_count):
        clean_array.append([])
        if len(raw_MO_data[val+3].split()) == 9: #want val+3 because first three lines are not the data we want
            for i in [0, 1, 2, 3]:
                clean_array[val].append(raw_MO_data[val+3].split()[i])
        elif len(raw_MO_data[val+3].split()) != 9:
            clean_array[val].append(raw_MO_data[val+3].split()[0])
            clean_array[val].append('blank')
            clean_array[val].append('blank')
            clean_array[val].append(raw_MO_data[val+3].split()[1])
        else:
            break
    
    return clean_array

def determine_indexing_lines(MO_array):
    '''
    Breaks up the inputted MO_array into groups based on modulus 5, i.e. the MO_array = [1,2, 3, 4, 5,6,7,8,9,10,11,12]
    would output lines = [[1,2,3,4,5],[6,7,8,9,10],[11,12]]
    The reason it does this is because that's how the Gaussian output file is set up
    that is, the MOs are listed as 5 columns, then it jumps down to the next five
    '''
    
    MO = 0
    next_MO = 1
    counter = 0
    lines = [[MO_array[MO]]]
    
    while next_MO <= len(MO_array)-1:
        while (MO_array[MO]-1)%5 < (MO_array[next_MO]-1)%5 and abs(MO_array[MO] - MO_array[next_MO]) < 5:
            lines[counter].append(MO_array[next_MO])
            if  next_MO < len(MO_array)-1:
                next_MO = next_MO + 1
            else:
                break
        counter = counter + 1
        if MO_array[next_MO] == MO_array[-1] and lines[-1][-1] == MO_array[next_MO]:
            break
        lines.append([])
        MO = next_MO
        next_MO = MO + 1
        lines[counter].append(MO_array[MO])
        
    return lines
    


def cleanup_raw_coefs(raw_MO_data, line_count, MO_array, clean_array, lines):
    ''' Cleans up the raw MO data lines and outputs an ordered MO/AO array (AOs on rows, MOs on columns)'''
    
    jump_index = line_count + 3
    clean_index = 0
    temp_array = []
    temp_array2 = []            
    big_counter = 0
    lines_max = len(lines)

    while big_counter < lines_max:
        if str(lines[big_counter][0]) in raw_MO_data[clean_index].split():
            for MO in lines[big_counter]:
                MO_col = locate_MO_index(MO,raw_MO_data[clean_index])
                for i in range(line_count):
                    temp_array.append([])
                    temp_array[i].append(raw_MO_data[clean_index+i+3].split()[-1:-6:-1])
                    temp_array2.append([])
                    temp_array2[i]=temp_array[i][0][::-1]
                    clean_array[i].append(temp_array2[i][MO_col])
                temp_array=[]
                temp_array2 =[]
            big_counter = big_counter + 1
        else:
            clean_index = clean_index + jump_index
    
    data_array = clean_array
    
    return data_array
    

def N2_CSPA_contributions(data_array, MO_array, line_count):
    '''performs C squared population analysis (see Ros and Schuit, Theoretica chimica acta, 1966, 4, 1, 1-12)
    If you're not working with Ag6-N2, you'll need to change the second for loop to reflect the atom order'''
    
    N7 = np.zeros(len(MO_array))
    N8 = np.zeros(len(MO_array))
    tots_sq = np.zeros(len(MO_array))
    N7_sq_rat = np.ones(len(MO_array))
    N8_sq_rat = np.ones(len(MO_array))
    
    for i in data_array:
        k = 0
        for j in i[4::]:
            tots_sq[k] = tots_sq[k] + abs(float(j))**2
            k = k + 1
        
    for i in range(len(data_array)):
        if str(9) in data_array[i][1]: # must change these if the N order moves!!!
            N7 = sq_sum_N(i,data_array,N7, line_count, 9)
        if str(10) in data_array[i][1]:
            N8 = sq_sum_N(i,data_array,N8, line_count, 10)
        
    for i in range(len(tots_sq)):
        N7_sq_rat[i] = 100.0*N7_sq_rat[i]*N7[i]/tots_sq[i]
        N8_sq_rat[i] = 100.0*N8_sq_rat[i]*N8[i]/tots_sq[i]
    
    pretty_N7 = [round(x,4) for x in N7_sq_rat]
    pretty_N8 = [round(x,4) for x in N8_sq_rat]
    
    return pretty_N7, pretty_N8
    

def sq_sum_N(i, data_array, N, line_count, atom_num):
    skipper = 0
    while data_array[i+skipper][1] == 'blank' or data_array[i+skipper][1] == str(atom_num):
        n = 0
        for j in data_array[i+skipper][4::]:
            N[n] = N[n] + abs(float(j)) ** 2
            n += 1
        skipper += 1
        if (i+skipper) == line_count:
            break
    
    return N
    
    
def test_AO_tots(data_array):
    '''use with all MOs to test that AOs sum to 1'''
    k = 0
    row_sum_occ = np.zeros(len(data_array))
    row_sum_tot = np.zeros(len(data_array))
    
    for i in data_array:
        for j in i[4:68]:
            row_sum_occ[k] = row_sum_occ[k] + (float(j)) ** 2
        for m in i[4::]:
            row_sum_tot[k] = row_sum_tot[k] + (float(m)) ** 2
        k += 1

def N2_contributions(data_array, MO_array, line_count):
    
    N7 = np.zeros(len(MO_array))
    N8 = np.zeros(len(MO_array))
    tots = np.zeros(len(MO_array))
    tots_sq = np.zeros(len(MO_array))
    N7_rat = np.ones(len(MO_array))
    N8_rat = np.ones(len(MO_array))

    
    for i in data_array:
        k = 0
        horiz = 0.0
        horiz_sq = 0.0
        for j in i[4::]:
            tots[k] = tots[k] + float(j)
            tots_sq[k] = tots[k] + abs(float(j))**2
            k = k + 1
            horiz = horiz + float(j)
            horiz_sq = horiz_sq + float(j) ** 2
        
    for i in range(len(data_array)):
        if str(9) in data_array[i][1]: # must change these if the N order moves!!!
            N7 = sum_N(i,data_array,N7, line_count, 9)
        if str(10) in data_array[i][1]:
            N8 = sum_N(i,data_array,N8, line_count, 10)
        

    for i in range(len(tots)):
        N7_rat[i] = N7_rat[i]*N7[i]/tots[i]
        N8_rat[i] = N8_rat[i]*N8[i]/tots[i]
    
def sum_N(i, data_array, N, line_count, atom_num):
    skipper = 0
    while data_array[i+skipper][1] == 'blank' or data_array[i+skipper][1] == str(atom_num):
        n = 0
        for j in data_array[i+skipper][4::]:
            N[n] = N[n] + float(j)
            n += 1
        skipper += 1
        if (i+skipper) == line_count:
            break
    
    return N
