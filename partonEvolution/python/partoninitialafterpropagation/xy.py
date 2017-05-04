#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np


fileobj = open("../ampt_0_0/ampt/ana/parton-initial-afterPropagation.dat")

multi = int(fileobj.readline().split()[2])
forwardxlist = []
forwardylist = []
backwardxlist = []
backwardylist = []
forwardNum = 0
backwardNum = 0

for i in range(multi):
    eachlineList = fileobj.readline().split()
    pz = float(eachlineList[3])
    if pz > 0:
        forwardxlist.append(float(eachlineList[5]))
        forwardylist.append(float(eachlineList[6]))
        forwardNum += 1
    else:
        backwardxlist.append(float(eachlineList[5]))
        backwardylist.append(float(eachlineList[6]))
        backwardNum += 1

print forwardNum, backwardNum
fig = plt.figure(figsize=(10,10), dpi=100)
fig.set_size_inches(6, 6)
ax1 = plt.scatter(forwardxlist, forwardylist, s=30, marker='o', c='r')
ax2 = plt.scatter(backwardxlist, backwardylist, s=30, marker='o', c='b')
plt.xlim(-8, 8)
plt.ylim(-8, 8)
plt.xlabel("x [fm]")
plt.ylabel("y [fm]")
plt.legend((ax1, ax2), ("pz > 0", "pz < 0"))
plt.show()
