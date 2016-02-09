FFLAGS = -g
#FFLAGS = -g -nbs -C -d2 -FR
#FFLAGS = -g -nbs -C 
exe    = gdmc

FC  = gfortran
#FC  = ifc -Vaxlib
#FC  = g77
INCL = 
MPIFLAG = 
LIBS = 
###############################################################
# MPI parameters: replace identical parameters if MPI needed  #
###############################################################
#FC  = pgf90 -Mfree
#INCL = -I/opt/scali/include 
#MPIFLAG = -DMPI
#LIBS    = -L/opt/scali/lib -lfmpi -lmpi
################## end of MPI section #########################

OBJ =  kmc.o io.o trans.o get.o cutstr.o

#################################################################################
#  targets
##################################################################################
# goal

$(exe):	$(OBJ) 
	$(FC) $(FFLAGS) -o $(exe) $(OBJ) $(LIBS)

# targets

kmc.o : kmc.f
	$(FC) $(FFLAGS) -c $*.f
io.o : io.f
	$(FC) $(FFLAGS) -c $*.f
trans.o: trans.f
	$(FC) $(FFLAGS) -c $*.f
get.o: get.f
	$(FC) $(FFLAGS) -c $*.f
cutstr.o: cutstr.f
	$(FC) $(FFLAGS) -c $*.f

clean:
	@rm -f $(OBJ) $(exe)

