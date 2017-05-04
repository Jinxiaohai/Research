#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Circle


fileobj = open("../../ampt/ana/npart-xy.dat")
xlist = []
ylist = []
eventHead = fileobj.readline()
multiProjectile = int(eventHead.split()[2])
multiTarget = int(eventHead.split()[3])

projectileX = []
projectileY = []
targetX = []
targetY = []
partPX = []
partPY = []
obsPx = []
obsPy = []
partTX = []
partTY = []
obsTx = []
obsTy = []

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
        if status > 0:
            partPX.append(x)
            partPY.append(y)
        else:
            obsPx.append(x)
            obsPy.append(y)
        
    if iteration < 0:
        targetX.append(x)
        targetY.append(y)
        if status > 0:
            partTX.append(x)
            partTY.append(y)
        else:
            obsTx.append(x)
            obsTy.append(y)


fig = plt.figure(figsize=(10,10), dpi=150)
fig.set_size_inches(6, 6)

plt.scatter(obsPx, obsPy, s=100, c='w', marker='o')
ax1 = plt.scatter(partPX, partPY, s=100, c='r', marker='o')
plt.scatter(obsTx, obsTy, s=100, c='w', marker='o')
ax2 = plt.scatter(partTX, partTY, s=100, c='g', marker='o')

plt.xlim(-10, 12)
plt.ylim(-10, 10)
plt.xlabel("x [fm]")
plt.ylabel("y [fm]")
plt.legend((ax1, ax2), ("Au", "Cu"))
plt.savefig("initial.png")
plt.show()
