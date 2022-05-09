#=============================================================================
#
#  Makefile for ipns (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 3-4-96
#
#=============================================================================
NPROC  = 4
NAME   = ipns 
DEBUG  = 
FFLAGS = -O -qrealsize=8 -qfree=f90 -c $(DEBUG)
OFLAGS = -O -qrealsize=8 $(DEBUG) -o $(NAME)
#LIB    = -lcomplib.sgimath
COMP   = xlf90
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  These objects depend on stuff.f90
#
MODS = consts.o  diff.o 
OBJS = ipns.o grad1.o grad2.o error.o
OBJECTS = rtflsp.o exit.o
#
#  End of objects
#
$(NAME): $(MODS) $(OBJS) $(OBJECTS)
	 $(COMP) $(OFLAGS) $(MODS) $(OBJS) $(OBJECTS) $(LIB)

clean:
	/bin/rm *.o *.kmo

$(OBJS): consts.o diff.o

$(MODS):

.f90.o:
	cp $*.f90 $*.F
	$(COMP) $(FFLAGS) $*.F 
	rm $*.F

.f.o:
	f77 -col120 $(FFLAGS) $*.f
