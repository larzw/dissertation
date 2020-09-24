EXE=out_iso1
FCFLAGS=-O2 -mcmodel=large -fopenmp

all:
	gfortran -llapack $(FCFLAGS) gauss_legendre.f90 inputs.f90 bonn.f90 ab.f90 x.f90 main.f90 -o $(EXE)

clean:
	rm -f *.mod

veryclean: clean
	rm -f *$(EXE)

