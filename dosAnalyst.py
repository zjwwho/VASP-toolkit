# -*- coding: utf-8 -*-
''' Module for VASP density of states data(DOSCAR file) analysis
    Author: https://github.com/zjwwho/VASP-toolkit
'''
from internalTools import DOS, goodsplit, writeDOS2file
from functools import reduce

def parseDOSCAR(DOSCAR_file, ldos_atoms):
    # function that parse DOSCAR file
    # ldos_atoms: (list type that contains atom No.)local-DOS of selected atoms
    fin = open(DOSCAR_file, 'r')
    data = fin.readlines()
    # obtain number of atoms
    atomnum = int(goodsplit(data[0])[0])
    # obtain energy range and NEDOS
    emax, emin, nedos = tuple(goodsplit(data[5])[:3])
    emax = float(emax)
    emin = float(emin)
    nedos = int(nedos)
    # is the calculation consider spin polarization
    ispin = True if len(goodsplit(data[6]))==5 else False
    del data[:5]
    
    # datablock include energy and dos
    datablock = []
    for i in range(atomnum+1):
        datablock.append(data[1:nedos+1])
        del data[:nedos+1]
    # total dos
    tdos = datablock[0]
    del datablock[0]
    tdos = [list(map(float, goodsplit(i))) for i in tdos]
    tdos = [[tdos[j][i] for j in range(len(tdos))] for i in range(len(tdos[0]))]
    tdos = DOS(tdos[0], tdos[1:], ispin, True)
    # partial dos
    pdos = [[list(map(float, goodsplit(i))) for i in datablock[j]] for j in range(len(datablock))]
    pdos = [[[pdos[i][k][j] for k in range(len(pdos[0]))] for j in range(len(pdos[0][0]))] for i in range(len(pdos))]
    pdos = [DOS(m[0], m[1:], ispin, False) for m in pdos]
    
    # write total dos
    tfout = open(DOSCAR_file+'.tdos.dat', 'w')
    writeDOS2file(tfout, tdos)
    
    # write local dos
    ldos = reduce(lambda x, y: x+y, [pdos[i] for i in ldos_atoms])
    lfout = open(DOSCAR_file+'.ldos.dat', 'w')
    writeDOS2file(lfout, ldos)
    
    print('parse DOSCAR finished!')