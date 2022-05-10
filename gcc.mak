#=============================================================================
#
#  Makefile for ipns (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 3-4-1996
#  Revised: 1-17-2020
#
#=============================================================================
NPROC  = 4
NAME   = ipns 
DEBUG  = 
OPT    = -O2
FFLAGS = -cpp -ffixed-line-length-120 -fdefault-real-8 -fdefault-double-8 \
-std=legacy $(DEBUG)
F90FLAGS = -cpp -fdefault-real-8 -fdefault-double-8 $(DEFINES) $(OPT) $(DEBUG)
OFLAGS = $(DEBUG) $(OPT) -o $(NAME)
LIB    = -L$(HOME)/local/OpenBLAS/lib -lopenblas
FC     = gfortran
F77    = gfortran
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  These objects depend on stuff.f90
#
MODS = consts.o diff.o 
OBJS = ipns.o grad1.o grad2.o error.o bslib1.o bslib2.o
OBJECTS = exit.o
#
# Optionally use Numerical-Recipes (commercial licensed)
#
ifdef USE_NR
  ifeq ($(LIBNR_DIR),)
    LIBNR_DIR = $(HOME)/git/NR-utilities
  endif
  LIB += -L$(LIBNR_DIR) -lnr
endif
#ifdef USE_NR
#  OBJECTS += nr_rtflsp.o
#endif
#
#  End of objects
#
$(NAME): $(MODS) $(OBJS) $(OBJECTS) 
	$(FC) $(OFLAGS) $(MODS) $(OBJS) $(OBJECTS) $(LIB)

clean:
	$(RM) -f *.mod *.o fort.* $(NAME)

.f90.o:
	$(FC) $(F90FLAGS) -c $*.f90 

.f.o:
	$(F77) $(FFLAGS) -c $*.f
