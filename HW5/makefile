CXX=g++-10
CC=gcc-10
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas -lfftw3 -lblas

CDEPS = nrutil.h auxiliary_functions.h

CPPDEPS = blas.h fields.hh fluid_2d.hh common.hh mgs_fem.hh tgmg.hh tgmg_base.hh tgmg_predict.hh tgmg_config.hh

OBJS1 = main_prob1.o nrutil.o auxiliary_functions.o fluid_2d.o common.o mgs_fem.o tgmg.o tgmg_base.o tgmg_predict.o tgmg_config.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

all: prob1_solve

prob1_solve: $(OBJS1)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob1: prob1_solve

clean:
	rm -f prob1_solve
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
