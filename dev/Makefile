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

OBJS= constants.o matrix-operations.o chebyshev.o model-canonical.o solver.o driver.o

instanton: %: $(OBJS)
	$(FC) $(FFLAGS) -o instanton $(OBJS) $(FLIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ $(FLIBS)

.PHONY : clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -f instanton
