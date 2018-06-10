# Default to GNU compiler, without warnings or profiling
FC = gfortran

WARN = T
PROF = F

ifeq ($(FC),gfortran)
	FFLAGS= -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp
ifeq ($(WARN),T)
	FFLAGS+= -Wall -fbounds-check -pedantic -fbacktrace
endif
ifeq ($(PROF),T)
	FFLAGS+= -g -pg
endif
endif

FLIBS= -lm -llapack -lblas

OBJS_UT = constants.o utils.o matrix-operations.o chebyshev.o

OBJS = $(OBJS_UT) model-canonical.o solver.o instanton.o driver.o
OBJS_OO = $(OBJS_UT) model-canonical.o solver.o instanton.o driver_instanton.o

instanton_oo: %: $(OBJS_OO)
	$(FC) $(FFLAGS) -o instanton $(OBJS_OO) $(FLIBS)

instanton: %: $(OBJS)
	$(FC) $(FFLAGS) -o instanton $(OBJS) $(FLIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ $(FLIBS)

.PHONY : clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -f instanton

cleandat:
	rm -f actions.dat instanton.dat solver-output.dat instanton_interp.dat instanton_.dat
