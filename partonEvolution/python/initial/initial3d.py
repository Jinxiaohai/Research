#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


fileobj = open("../ampt_0_0/ampt/ana/npart-xy.dat")
xlist = []
ylist = []
eventHead = fileobj.readline()
multiProjectile = int(eventHead.split()[2])
multiTarget = int(eventHead.split()[3])

projectileX = []
projectileY = []
projectileZ = []
targetX = []
targetY = []
targetZ = []

partPX = []
partPY = []
partPZ = []
obsPx = []
obsPy = []
obsPz = []

partTX = []
partTY = []
partTZ = []
obsTx = []
obsTy = []
obsTz = []

for inum in range(multiProjectile+multiTarget):
    eachline = fileobj.readline().split()
    x = float(eachline[0])
    y = float(eachline[1])
    iteration = int(eachline[2])
    status = int(eachline[3])
    z = float(eachline[4])
    if iteration > 0:
        projectileX.append(x)
        projectileY.append(y)
        projectileZ.append(z)
        if status > 0:
            partPX.append(x)
            partPY.append(y)
            partPZ.append(z)
        else:
            obsPx.append(x)
            obsPy.append(y)
            obsPz.append(z)
        
    if iteration < 0:
        targetX.append(x)
        targetY.append(y)
        targetZ.append(z)
        if status > 0:
            partTX.append(x)
            partTY.append(y)
            partTZ.append(z)
        else:
            obsTx.append(x)
            obsTy.append(y)
            obsTz.append(z)


fig = plt.figure(projection='3d')
fig.set_size_inches(14, 14)

plt.scatter(obsPx, obsPy, obsPz, s=60, c='w', marker='o')
plt.scatter(partPX, partPY, partPZ, s=60, c='r', marker='o')
plt.scatter(obsTx, obsTy, obsTz, s=60, c='w', marker='o')
plt.scatter(partTX, partTY, partTZ, s=60, c='g', marker='o')

plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.zlim(-10, 10)
plt.show()
