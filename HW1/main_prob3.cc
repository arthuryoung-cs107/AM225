#include <cstdio>
#include <cstdlib>
#include <stdint.h>

#include "HW1_aux.hh"
#include "omp.h"
#include "gsl/gsl_rng.h"

extern "C" {
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

int main()
{
  long check_prime = 100;
  primes * case1 = find_primes(check_prime);

  return 0;
}
