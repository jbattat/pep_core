FORT = gfortran

FFLAGS = -fno-automatic -w -frecord-marker=8
B = ../bigtest
G = .
P = ../pep
U = ../peputil

new:
	$(FORT) $(FFLAGS) -c *.f

# include dependencies:
aplot.o: $U/plttyp.inc $G/pltopnfo.inc
symbol.o: $U/plttyp.inc
