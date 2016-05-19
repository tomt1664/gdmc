#!/usr/bin/python

# histogram.py
# T.Trevethan
# 2015
#
# Script to calculate a histogram of mean squared displacement (MSD)
# by averaging over the trajectory of a diffusing defect or atom
#
# Defect displacement and times are taken from the dist.dat file
# created by the displace.py script
#
# Two inputs are required: the histogram size, or maximum time: set with the -m flag
#                          and the bin size: set with the -b flag 

import sys
import math 
import argparse
import numpy as np

# Parse command line input
parser = argparse.ArgumentParser()
parser.add_argument('-m', type = float, help = 'histogram size [float]')
parser.add_argument('-b', type = float, help = 'bin size [float]')
variables = parser.parse_args()
maxt = variables.m
bint = variables.b

if not maxt:
    print "Error: -m flag not set"
    sys.exit(1)

if not bint:
    print "Error: -b flag not set"
    sys.exit(1)

# Open input, time and output files
input  = open("dist.dat", 'r')
hxout = open("hist_x.dat",'w') 
hyout = open("hist_y.dat",'w')

#read input file and store values in array
displ = []

for line in input:
    fdata = line.split()
    data = (float(fdata[0]),float(fdata[1]),float(fdata[2]),float(fdata[3]),float(fdata[4]))
    displ.append(data)

print maxt
print bint

print "No. of entries: " + str(len(displ))

# number or runs to average over
nruns = int(displ[-1][0]/maxt)
print "No. samples: " + str(nruns)

# number of bins
nbins = int(maxt/bint) + 1
print "No. bins: " + str(nbins)

tadd = 0.0
irpt = 0
daad1 = 0.0
daad2 = 0.0

# displacement bin arrays
htx = np.zeros((nbins,nruns))
hty = np.zeros((nbins,nruns))

# squared displacement bin arrays
htxs = np.zeros((nbins,nruns))
htys = np.zeros((nbins,nruns))

for i in range(nruns):
    nbin = 0
    numbin = 0

    for j in range(len(displ)):
        time = displ[j+irpt][0] - tadd

        if time < nbin*bint:
            htx[nbin][i] = htx[nbin][i] + displ[j+irpt][1] - daad1
            hty[nbin][i] = hty[nbin][i] + displ[j+irpt][2] - daad2

            htxs[nbin][i] = htxs[nbin][i] + (displ[j+irpt][1] - daad1)**2
            htys[nbin][i] = htys[nbin][i] + (displ[j+irpt][2] - daad2)**2    
            
            numbin += 1
            
        else:
            if numbin == 0:
                htx[nbin][i] = 0.0
                hty[nbin][i] = 0.0
                htxs[nbin][i] = 0.0
                htys[nbin][i] = 0.0
            else:
                htx[nbin][i] = htx[nbin][i]/(numbin*1.0)
                hty[nbin][i] = hty[nbin][i]/(numbin*1.0)
                htxs[nbin][i] = htxs[nbin][i]/(numbin*1.0)
                htys[nbin][i] = htys[nbin][i]/(numbin*1.0)

            numbin = 0
            nbin += 1
        
        if nbin*bint >= maxt: 
            jend = j
            break
    
    tadd = displ[j+irpt][0]
    daad1 = displ[j+irpt][1]
    daad2 = displ[j+irpt][2]
    irpt += jend

# add together the binned runs
histx = np.zeros(nbins)
histy = np.zeros(nbins)
histxsq = np.zeros(nbins)
histysq = np.zeros(nbins)

for ii in range(nbins-1):
    for jj in range(nruns):
        histx[ii] += htx[ii,jj]/(nruns*1.0)
        histy[ii] += hty[ii,jj]/(nruns*1.0)
        histxsq[ii] += htxs[ii,jj]/(nruns*1.0)
        histysq[ii] += htys[ii,jj]/(nruns*1.0)

    outputlinex = str(ii*bint)+" "+str(histx[ii])+" "+str(histxsq[ii])+"\n"
    hxout.write(outputlinex)
    outputliney = str(ii*bint)+" "+str(histy[ii])+" "+str(histysq[ii])+"\n"
    hyout.write(outputliney)

print "complete"

