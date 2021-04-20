CXX=g++-10
CC=gcc-10
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas -lfftw3 -lblas

CDEPS = nrutil.h auxiliary_functions.h knuth_lcg.h

CPPDEPS = poisson_fft_AY.hh blas.h Ritz_Galerk.hh conj_grad.hh quadrat.hh

OBJS2 = main_prob2.o poisson_fft_AY.o nrutil.o auxiliary_functions.o knuth_lcg.o Ritz_Galerk.o conj_grad.o quadrat.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

all: prob2_solve

prob2_solve: $(OBJS2)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob2: prob2_solve

clean:
	rm -f prob2_solve
	rm -f *.o
clean_dat:
	rm -f *.dat
	rm -f *.aysml
	rm -f *.aydat
clean_all: clean clean_dat
