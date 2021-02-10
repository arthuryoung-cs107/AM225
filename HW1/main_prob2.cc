#include <cstdio>
#include <cstdlib>
#include <stdint.h>

#include "HW1_aux.hh"
#include "omp.h"

extern "C" {
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

int main()
{
  int i, j;
  int m = 80;
  int n = 40;
  int generations = 150;
  uint64_t seed_start = 88888888;

  int ** cell = imatrix(0, m-1, 0, n-1);
  int ** cell_old = imatrix(0, m-1, 0, n-1);

  cell_automaton(m, n, 1);

}
