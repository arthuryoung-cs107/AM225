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
  long check_prime = 1e5;
  Primes * case1 = find_primes(check_prime);

  // long power_1 = 8258993;
  // long m_base = 40;

  long power_1 = 7;
  long m_base = 2;

  Mersenne * test1 = Mersenne_expand(power_1, m_base);
  long remain = general_divide(test1, 2);

  
  long num_primes = Count_Primes();

  free_Primes(case1);
  free_Mersenne(test1);
  return 0;
}
