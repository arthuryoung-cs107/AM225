CXX=g++-10
CC=gcc-10
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas -lfftw3 -lblas

CDEPS = nrutil.h auxiliary_functions.h knuth_lcg.h

CPPDEPS = blas.h AYmat.hh ADM.hh

OBJS = main.o nrutil.o auxiliary_functions.o knuth_lcg.o AYmat.o ADM.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

all: RPCA_test

RPCA_test: $(OBJS)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

RPCA: RPCA_test

clean:
	rm -f RPCA_test
	rm -f *.o
clean_dat:
	rm -f *.dat
	rm -f *.aysml
	rm -f *.aydat
clean_datdir:
	rm -f ./dat_dir/*.dat
	rm -f ./dat_dir/*.aysml
	rm -f ./dat_dir/*.aydat

clean_all: clean clean_dat
