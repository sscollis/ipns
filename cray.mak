#=============================================================================
#
#  Makefile for ipns (Cray)
#
#  Author:  Scott Collis
#
#  Revised: 3-4-96
#
#=============================================================================
NPROC  = 4
NAME   = ipns 
DEBUG  = -O1
FFLAGS = -c $(DEBUG)
OFLAGS = $(DEBUG) -o $(NAME)
LIB    = -limsl
COMP   = f90
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  These objects depend on stuff.f90
#
MODS = consts.o  diff.o 
OBJS = ipns.o grad1.o grad2.o error.o
OBJECTS = rtflsp.o
#
#  End of objects
#
$(NAME): $(MODS) $(OBJS) $(OBJECTS)
	 $(COMP) $(OFLAGS) $(MODS) $(OBJS) $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o

$(OBJS): consts.o diff.o
	$(COMP) $(FFLAGS) -p consts.o -p diff.o $*.f90

$(MODS):
	$(COMP) $(FFLAGS) $*.f90

.f90.o:
	 $(COMP) $(FFLAGS) $*.f90 

.f.o:
	 cf77 -N80 $(FFLAGS) $*.f
