CXX=g++-10
CC=gcc-10
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas


CDEPS = nrutil.h auxiliary_functions.h

CPPDEPS = Cash_Karp.hh HW2_aux.hh

OBJS2 = main_prob2.o nrutil.o auxiliary_functions.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

prob2_solve: $(OBJS2)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob2: prob2_solve

all: prob2_solve

clean:
	rm -f prob2_solve
	rm -f *.o
clean_dat:
	rm -f *.dat
	rm -f *.aysml
	rm -f *.aydat
clean_all: clean clean_dat
