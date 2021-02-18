#ifndef HW1_AUX_HH
#define HW1_HEADER_HH

#include <cstdio>
#include <stdint.h>
#include <math.h>
#include "gsl/gsl_rng.h"

struct Prime_struct
{
  long * primes;
  long max;
  long len;
} ;
typedef struct Prime_struct Primes;

struct Mersenne_Struct
{
  long * coeffs;
  long k_max;
  long m;
  long base;
} ;
typedef struct Mersenne_Struct Mersenne;

double random_uni(double low, double high, uint64_t * carry);
uint64_t casino_game(gsl_rng * T);
void cell_automation(int m, int n, int gen_max, double * out_vec);
Primes * find_primes(long check_prime);
Mersenne * Mersenne_expand(long pow, long m_base);
long general_divide(Mersenne * num_in, long divide);
void free_Primes(Primes * prime_structure);
void free_Mersenne(Mersenne * mers_struct);

#endif
