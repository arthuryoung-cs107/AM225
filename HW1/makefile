CXX=g++-10
CC=gcc-10
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas


CDEPS = nrutil.h auxiliary_functions.h knuth_lcg.h

CPPDEPS = HW1_aux.hh

OBJS1 = main_prob1.o nrutil.o auxiliary_functions.o knuth_lcg.o HW1_aux.o
OBJS2 = main_prob2.o nrutil.o auxiliary_functions.o knuth_lcg.o HW1_aux.o
OBJS2 = main_prob3.o nrutil.o auxiliary_functions.o knuth_lcg.o HW1_aux.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

prob1_solve: $(OBJS1)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob2_solve: $(OBJS2)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob3_solve: $(OBJS3)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob1: prob1_solve

prob2: prob2_solve

all: prob1_solve prob2_solve prob3_solve

clean:
	rm -f prob1_solve
	rm -f prob2_solve
	rm -f prob3_solve
	rm -f *.o
clean_dat:
	rm -f *.dat
clean_all: clean clean_dat
