# -*- coding: utf-8 -*-
"""
Created on Fri Sep 02 11:52:34 2016

@author: CHILDRESSLAB
"""

import numpy as np
import csv
from ConexXYZMove import *

f = open('MagnetPath.csv', 'r')
MagPath = csv.reader(f)

alpha = np.pi/2
x0 = 14
y0 = 13
z0 = -4

def Rtrans(R, a, x0, y0, z0):
    """
    Transforms input vector R to new vector R2 by first rotating by a about
    z axis, and then translating by x0, y0, z0.
    """
    
    #define rotation matrix about z
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Rrot = np.dot(Rz,np.array(R))
    Rtrans = np.array([x0,y0,z0])
    
    R2 = Rrot + Rtrans
    
    return R2
    
    
for row in MagPath:
    R = Rtrans(map(float,row), alpha, x0, y0, z0)
    SetPosition(R)
    WaitForMove()
    print GetPosition()