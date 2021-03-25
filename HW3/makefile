CXX=g++-10
CC=gcc-10
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas

CDEPS = nrutil.h auxiliary_functions.h

CPPDEPS = HW3_aux.hh

OBJS5 = main_prob6.o nrutil.o auxiliary_functions.o HW3_aux.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

all: prob5_solve

prob5_solve: $(OBJS5)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob5: prob5_solve

clean:
	rm -f prob5_solve
	rm -f *.o
clean_dat:
	rm -f *.dat
	rm -f *.aysml
	rm -f *.aydat
clean_all: clean clean_dat
