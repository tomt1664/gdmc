# gdmc
Kinetic Monte Carlo code for simulating point defect dynamics and aggregation in graphitic systems

GraphDMC is an on-lattice Kinetic Monte Carlo/Dynamic Monte Carlo (KMC) code for simulating real-time point defect evolution in the graphite (and implicitly graphene) structure, written in Fortran. Point defects in this context are vacancies (atoms missing from the crystal lattice) and interstitials (extra atoms between the planes), not topological defects resulting from bond rotations such as Stone-Wales defects. 

The code employs the standard KMC or 'n-fold way' method which is an efficient method for simulating the dynamical evolution, and hence enabling detailed sampling, of complex processes consisting many individual transitions. In essence, the algorithm is used to decide which out of a number of different possible transitions, each with different rates, is followed to leave a particular state, as well as determining the time spent residing in that state. 

This repository contains the source code files, makefile, and an examples directory containing example input and potential energy surface files. 

The code has been tested with the gfortran compiler, which must be available from the command line. The compiler can be changed by editing the Makefile. It is not guaranteed to work for other compilers, but it shouldn't be a problem.  To compile the code, simply type:

> $ make

This creates an executable in the same directory called gdmc 

The program is then run in the command line. It requires no standard input, but writes to standard output, which can be piped into an output file:

> $ ./gdmc > outputfile

The code requires four formatted files in the directory in which it is run: gkin.dat, ens.dat, iens.dat and system.in

In addition to the standard output, it creates three output files by default: tens.out, system.out and system.xyz. Depending on the input options it can also create the files: vacs.xyz, ints.xyz and link.xyz, which give vacancy, interstitial and interlayer di-vacancy coordinates respectively. 

For the full documentation, see theattached wiki. 

Examples of usage: <a href="http://pubs.rsc.org/en/content/articlehtml/2014/nr/c3nr06222h">Strain directed diffusion</a>

##Utilities

The pytools directory contains a collection of Python scripts to perform analysis and processing of the code output. 

> displace.py

This script calculates the displacement of a single defect as a function of simulation time. It takes information from the system.xyz trajectory file and the tens.out file. For correct usage set 
