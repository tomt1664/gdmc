#!/usr/bin/python

# displace.py
# T.Trevethan
# 2015
#
# Script to calculate the displacement of a single defect during a gdmc simulation
#
# Defect coordinates are taken from the system.xyz trajectory file and simulation time from the tens.out file
#
# Periodic cell dimensions (system matrix) are supplied with -x, -y, and -l arguments
#
# The time and x,y,z components and total displacement are written to the dist.dat file 

import sys
import math 
import argparse

#set open boundary conditions if not passed
xdim = 10000
ydim = 10000
zdim = 10000

# Parse command line input
parser = argparse.ArgumentParser()
parser.add_argument('-x', type = int, help = 'x periodic dimenstion [integer]')
parser.add_argument('-y', type = int, help = 'y periodic dimenstion [integer]')
parser.add_argument('-l', type = int, help = 'number of layers [integer]')
variables = parser.parse_args()
xdim = 1.22*variables.x
ydim = 0.71*variables.y
zdim = 3.35*variables.l

# Open input, time and output files
input  = open("system.xyz", 'r')
time = open("tens.out",'r') 
output = open("dist.dat",'w')

# Read first coordinate
numat = input.readline()
nat = int(numat)

# trajectory file should have only one defect
if nat != 1:
    print "Error: More than one defect in the coordinate file"
    sys.exit(1)

blank = input.readline()
atline = input.readline()
atdata = atline.split()
icoord = [float(atdata[1]), float(atdata[2]), float(atdata[3])]

print "Initial position: " + str(icoord[0]) + " " + str(icoord[1]) + " " + str(icoord[2])
outputline = "0.0 0.0 0.0 0.0\n"
output.write(outputline)

#keep track of boundary crossings
xadd = 0
yadd = 0
zadd = 0
coord = icoord

line = input.readline()
#loop over all lines in input file
while line != "":
    nat = int(line)
    if nat == 1:
        #keep old coordinates
        ocoord = coord
        blank = input.readline()    
        atline = input.readline()
        atdata = atline.split()
        coord = [float(atdata[1]), float(atdata[2]), float(atdata[3])]

        #track boundary crossings
        if coord[0] < xdim*0.3 and ocoord[0] > xdim*0.7:
            xadd = xadd - xdim
        if ocoord[0] < xdim*0.3 and coord[0] > xdim*0.7:
            xadd = xadd + xdim
        if coord[1] < ydim*0.3 and ocoord[1] > ydim*0.7:
            yadd = yadd - ydim
        if ocoord[1] < ydim*0.3 and coord[1] > ydim*0.7:
            yadd = yadd + ydim
        if coord[2] < zdim*0.3 and ocoord[2] > zdim*0.7:
            zadd = zadd - zdim
        if ocoord[2] < zdim*0.3 and coord[2] > zdim*0.7:
            zadd = zadd + zdim

        #calculate displacement
        xdis = coord[0] - icoord[0] - xadd
        ydis = coord[1] - icoord[1] - yadd
        zdis = coord[2] - icoord[2] - zadd
        mdis = math.sqrt(xdis*xdis + ydis*ydis)
        tline = time.readline()
        tdata = tline.split()
        simtime = float(tdata[0])
        #write the components and magnitude to the output file
        outputline = str(simtime)+" "+str(xdis)+" "+str(ydis)+" "+str(zdis)+" "+str(mdis)+"\n"
        output.write(outputline)

    line = input.readline()

print "Final position: " + str(coord[0]) + " " + str(coord[1]) + " " + str(coord[2])
