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
  uint64_t check_prime = 1e5;
  uint64_t * prime_vec = find_primes(check_prime)

  return 0;
}
