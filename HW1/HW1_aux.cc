#include <cstdio>

#include "HW1_aux.hh"
#include "omp.h"

extern "C" {
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}


uint64_t casino_game(uint64_t * seed)
{
  double end_cond = 1;
  double cond = 0;
  uint64_t count = 0;

  while(cond < end_cond)
  {
    cond = cond + uniform_random(0, 1, seed);
    count++;
  }

  return count;
}
