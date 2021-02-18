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
  long check_prime = 2e5;
  Primes * case1 = find_primes(check_prime);

  long power_1 = 7;
  long m_base = 2;

  Mersenne * test1 = Mersenne_expand(power_1, m_base);
  long remain = general_divide(test1, 2);

  long power_2 = 8258993;
  long m_base2 = 40;
  Mersenne * M = Mersenne_expand(power_2, m_base2);
  long num_primes = Count_Primes(case1, M);
  printf("number of dividing primes: %ld \n", num_primes);

  long check_prime2 = 1e6;
  Primes * case2 = find_primes(check_prime2);

  long power_3 = 8258992;
  long m_base3 = 40;
  Mersenne * M2 = Mersenne_expand(power_3, m_base3);
  long num_primes2 = Count_Primes(case2, M2);
  printf("number of dividing primes: %ld \n", num_primes2);


  free_Primes(case1);
  free_Mersenne(test1);
  return 0;
}
