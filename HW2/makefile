CXX=g++-10
CC=gcc-10
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

LDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas

CDEPS = nrutil.h auxiliary_functions.h

CPPDEPS = Cash_Karp.hh Geng.hh HW2_aux.hh

OBJS2 = main_prob2.o nrutil.o auxiliary_functions.o Cash_Karp.o

OBJS6 = main_prob6.o nrutil.o auxiliary_functions.o HW2_aux.o Geng.o

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)

%.o: %.cc $(CPPDEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LDIR)

all: prob2_solve prob6_solve

prob2_solve: $(OBJS2)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob6_solve: $(OBJS6)
	$(CXX) -o $@ $^ $(CFLAGS) $(LINK)

prob2: prob2_solve

prob6: prob6_solve

clean:
	rm -f prob2_solve
	rm -f prob6_solve
	rm -f *.o
clean_dat:
	rm -f *.dat
	rm -f *.aysml
	rm -f *.aydat
clean_prob6_datdir:
	rm -f ./dat_dir/prob6_*

clean_all: clean clean_dat
