# -*- coding: utf-8 -*-
''' Module that contains internal tools for other modules
    Author: https://github.com/zjwwho/VASP-toolkit
'''
from math import sqrt

def goodsplit(line):
    # split string by space
    return [item for item in line.strip().split(' ') if item != '']

def ObliP2CartP(obliP, obliCS):
    # convert point in Oblique coordinate system to point in Cartesian coordinate system
    # obliP: point in Oblique coordinate system (Relative coordinates)
    # obliCS: vector of Oblique coordinate system
    obliCellPara = [sqrt(pow(obliCS[i][0],2)+pow(obliCS[i][1],2)+pow(obliCS[i][2],2)) for i in range(3)]
    cartP = [0.0, 0.0, 0.0]
    cartP[0] = (obliP[0]*(obliCS[0][0]/obliCellPara[0])+obliP[1]*(obliCS[1][0]/obliCellPara[1])+obliP[2]*(obliCS[2][0]/obliCellPara[2]))*obliCellPara[0]
    cartP[1] = (obliP[0]*(obliCS[0][1]/obliCellPara[0])+obliP[1]*(obliCS[1][1]/obliCellPara[1])+obliP[2]*(obliCS[2][1]/obliCellPara[2]))*obliCellPara[1]
    cartP[2] = (obliP[0]*(obliCS[0][2]/obliCellPara[0])+obliP[1]*(obliCS[1][2]/obliCellPara[1])+obliP[2]*(obliCS[2][2]/obliCellPara[2]))*obliCellPara[2]
    return cartP
    
def position2distance(position_list, cs=[[1,0,0],[0,1,0],[0,0,1]]):
    # convert position of points to distance
    dis = []
    pre = list(map(float, position_list[0]))
    for i in position_list:
        i = ObliP2CartP(list(map(float, i)), cs)
        dis.append(sqrt(pow(i[0]-pre[0], 2)+pow(i[1]-pre[1], 2)+pow(i[2]-pre[2], 2)))
        pre = i
    absdis = [dis[0]]
    for i in range(len(dis)-1):
        absdis.append(absdis[i]+dis[i+1])
    return absdis
    
class DOS():
    # class for density of states storage and add operation
    def __init__(self, energy, dos, ispin, istdos):
        self.energy = energy
        if istdos:
            if ispin:
                self.dos = dos[:-2]
            else:
                self.dos = dos[:-1]
        else:
            self.dos = dos
        self.ispin = ispin
        self.istdos = istdos
    # add
    def __add__(self, other):
        return DOS(self.energy, [[self.dos[m][n]+other.dos[m][n] for n in range(len(self.energy))] for m in range(len(self.dos))], self.ispin, self.istdos)
    
def writeDOS2file(fstream, dos):
    # write DOS type to file
    for i in range(len(dos.energy)):
        fstream.write('{:>10f}'.format(dos.energy[i]))
        for j in range(len(dos.dos)):
            fstream.write('  {:>10f}'.format(dos.dos[j][i]))
        fstream.write('\n')
    print('write DOS to file finished!')
    
def issamepoint(p1, p2, threshold=0.0001):
    # verify p1 and p2 are the same point
    flag = 0
    for i in range(len(p1)):
        if abs(p1[i]-p2[i]) <= threshold: flag += 1
    if flag == len(p1):
        return True
    else:
        return False
        
class POSCAR():
    # class contain crystal structure information(POSCAR/CONTCAR)
    def __init__(self, poscar):
        fin = open(poscar, 'r')
        data = [goodsplit(i) for i in fin.readlines()]
        # scaling factor
        sf = float(data[1][0])
        # lattice vector(multiply by scaling factor)
        self.lv = [list(map(lambda x: float(x)*sf, i)) for i in data[2:5]]
        # atom types
        self.at = data[5]
        # atom numbers
        self.an = list(map(int, data[6]))
        # is selective dynamics?
        self.sd = True if data[7][0].upper()=='S' else False
        if self.sd: del data[7]
        # cartesian or direct lattice
        self.cdl = 'Cartesian' if data[7][0].upper()=='C' or data[7][0].upper()=='K' else 'Direct'
        del data[:8]
        # atom positions
        self.ap = [list(map(float, i[:3])) for i in data[:sum(self.an)]]

    def getAtomNobyPosition(self, atomPosition):
        # get atom No. by atom position
        for i, p in enumerate(self.ap):
            if issamepoint(atomPosition, p):
                return i