# Load the common configuration file
include ../config.mk
CC=gcc-10
iflags=-I../tgmg
lflags=-L.

objs=common.o fluid_2d.o mgs_fem.o
src=$(patsubst %.o,%.cc,$(objs))
execs=prob1_solve

all:
	$(MAKE) -C ../tgmg
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

libf2d.a: $(objs)
	rm -f libf2d.a
	ar rs libf2d.a $^

c_lib.a: $(OBJS1)
	$(cxx) -o $@ $^

# fluid_test: fluid_test.cc libf2d.a
# 	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

%.o: %.c $(CDEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDIR)


prob1_solve: main_prob1.cc libf2d.a c_lib.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d

clean:
	rm -f $(execs) $(objs) libf2d.a

clean-all:
	$(MAKE) -C ../tgmg clean
	rm -f $(execs) $(objs) libf2d.a

.PHONY: clean all executables depend
