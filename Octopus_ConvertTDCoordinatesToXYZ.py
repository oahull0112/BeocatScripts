# -*- coding: utf-8 -*-
"""
@author: Olivia Hull
"""

def runOctParse(filename, tA, nA):
    '''
    Takes an Octopus coordinate file and outputs in .xyz Format
    inputs:
        filename: name of coordinate file (octopus names it "coordinates" by default)
        tA: an array of the atom types in order that they appear in the output file
        nA: the number of atoms of a type in order that they appear in the output file
            e.g if in the Octopus input file you specified:
                Ag | x1 | y1 | z1
                Ag | x2 | y2 | z2
                C  | x3 | y3 | z3
                H  | x4 | y4 | z4
                then tA = ["Ag", "C", "H"] and nA = [2, 1, 1]
            if you have something like:
                Ag | x1 | y1 | z1
                Ag | x2 | y2 | z2
                C  | x3 | y3 | z3
                H  | x4 | y4 | z4
                Ag | x5 | y5 | z5
                Ag | x6 | y6 | z6
                C  | x7 | y7 | z7
                then tA = ["Ag", "C", "H", "Ag", "C"]
                and  nA = [2, 1, 1, 2, 1]
            to troubleshoot that the program is getting it right, insert a print line
            on the variable anList to ensure it matches with your Octopus input file.
    Outputs:
        The program writes a file formatted as "inputfilename_xyzFormat.xyz"
        where inputfilename is the name you assigned to the variable filename
            
    '''
    nAtoms = sum(nA)
    outName = filename + '_xyzFormat.xyz'
    anList = getAtomNumberList(tA, nA)
    coords = getCoordinates(filename, nAtoms)
    makeOutputFile(coords, nAtoms, outName, anList)
    

def getAtomNumberList(tA, nA):
    anList = []
    c = 0
    for i in nA:
        for j in range(i):
            anList.append(tA[c])
        c = c+1
    return anList    


def getCoordinates(filename, nAtoms):
    with open(filename) as f:
        cList = f.readlines()
    
    cList = [x.strip() for x in cList[5::]] # clean up clist
    
    coords = []
    
    for i in range(len(cList)):
        coords.append(cList[i].split()[2:(nAtoms*3)+2]) # the [2::] is ignoring the first two entries, the iteration number and the time

    return coords 
#nAtoms = int((len(coords[0]))/3)

def makeOutputFile(coords, nAtoms, outName, anList):
    tList = []
    
    for j in range(len(coords)):
        tList.append(str(nAtoms))
        tList.append(' ')
        for i in range(nAtoms):
            i = i*3
            tList.append(coords[j][0+i:3+i])
    
    outF = open(outName, "w")
    
    for l in range(len(tList)):
        if l%(nAtoms+2) == 0:
            outF.write(tList[l])
            outF.write("\n")
        elif l%(nAtoms+2) == 1:
            outF.write("\n")
        else:
            j = l%(nAtoms + 2)
            q = j - 2
            aux_line = str(anList[q] + " " + str(tList[l][0]) + " " + str(tList[l][1]) + " " + str(tList[l][2]))
            outF.write(aux_line)
            outF.write("\n")
    
    outF.close()


