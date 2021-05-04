CXX=g++-10
CC=gcc-10
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas -lfftw3 -lblas

CDEPS = nrutil.h auxiliary_functions.h knuth_lcg.h

CPPDEPS = blas.h Ritz_Galerk_sphere.hh conj_grad.hh quadrat.hh AYmat.hh cubic_1d_fe.hh cubic_1d_alt_C2.hh cubic_1d_alt_C1.hh

OBJS2 = main_prob2.o nrutil.o auxiliary_functions.o knuth_lcg.o Ritz_Galerk_sphere.o conj_grad.o quadrat.o AYmat.o

OBJS3 = main_prob3.o nrutil.o auxiliary_functions.o cubic_1d_fe.o conj_grad.o quadrat.o cubic_1d_alt_C2.o cubic_1d_alt_C1.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

all: prob2_solve prob3_solve

prob2_solve: $(OBJS2)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob3_solve: $(OBJS3)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob2: prob2_solve

prob3: prob3_solve

clean:
	rm -f prob2_solve
	rm -f prob3_solve
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
