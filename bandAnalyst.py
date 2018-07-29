# -*- coding: utf-8 -*-
''' Module for VASP band structure data(EIGENVAL file) analysis
    Author: https://github.com/zjwwho/VASP-toolkit
'''
from internalTools import goodsplit, position2distance

def parseEIGENVAL(EIGENVAL_file, rl):
    # function that parse EIGENVAL file
    # rl: reciprocal lattice --- needed for kpoints distance calculation
    fin = open(EIGENVAL_file, 'r')
    data = fin.readlines()
    # obtain number of atoms
    atomnum = int(goodsplit(data[0])[0])
    # obtain number of kpoints and bands
    knum, bnum = tuple(map(int, goodsplit(data[5])[1:3]))
    # is the calculation consider spin polarization
    ispin = True if len(goodsplit(data[8]))==3 else False
    # cell volume
    cellV = float(goodsplit(data[1])[0])*atomnum
    # cell parameters
    cellP = tuple(map(lambda x: float(x)*1.0E10, goodsplit(data[1])[1:4]))
    del data[:6]

    # kblocks include kpoint position and bands
    kblocks = []
    for i in range(knum):
        kblocks.append(data[1:bnum+2])
        del data[:bnum+2]
    
    kpoints = [goodsplit(kblocks[i][0])[:3] for i in range(knum)]
    kpdis = position2distance(kpoints, rl)
    
    fout = open(EIGENVAL_file+'.dat', 'w')
    if ispin:
        # write spin-up and spin-down band
        for m in range(bnum):
            for n in range(knum):
                fout.write('{:>10f}  {:>10f}  {:>10f}\n'.format(kpdis[n], float(goodsplit(kblocks[n][m+1])[1]), float(goodsplit(kblocks[n][m+1])[2])))
            fout.write('\n')
    else:
        # write spin-up band
        for m in range(bnum):
            for n in range(knum):
                fout.write('{:>10f}  {:>10f}\n'.format(kpdis[n], float(goodsplit(kblocks[n][m+1])[1])))
            fout.write('\n')
    print('parse EIGENVAL file finished!')