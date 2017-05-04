#!/usr/bin/env python
# _*_ coding: utf-8 _*_

"""document of the module"""


import os
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np


fileobj = open("../../ampt/ana/after_parton.dat")

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
fig = plt.figure(figsize=(9.5,10), dpi=100)
plt.scatter(forwardxlist, forwardylist, s=40, marker='o', c='r')
plt.scatter(backwardxlist, backwardylist, s=40, marker='o', c='b')
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.xlabel("x [fm]")
plt.ylabel("y [fm]")
plt.savefig("zpc.png")
plt.show()
