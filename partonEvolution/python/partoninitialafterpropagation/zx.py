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
forwardzlist = []
backwardxlist = []
backwardzlist = []
forwardNum = 0
backwardNum = 0

for i in range(multi):
    eachlineList = fileobj.readline().split()
    pz = float(eachlineList[3])
    if pz > 0:
        forwardxlist.append(float(eachlineList[5]))
        forwardzlist.append(float(eachlineList[7]))
        forwardNum += 1
    else:
        backwardxlist.append(float(eachlineList[5]))
        backwardzlist.append(float(eachlineList[7]))
        backwardNum += 1

print forwardNum, backwardNum
fig = plt.figure(figsize=(9.5,10), dpi=100)
plt.scatter(forwardzlist, forwardxlist, s=60, marker='o', c='r')
plt.scatter(backwardzlist, backwardxlist, s=60, marker='o', c='b')
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.xlabel("z [fm]")
plt.ylabel("x [fm]")
plt.show()
