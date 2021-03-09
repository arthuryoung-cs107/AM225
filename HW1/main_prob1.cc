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
  uint64_t i, j, k, count, carry, thread_num;
  uint64_t trials = 1e9;
  uint64_t jump = 10;
  double mean_count, t0, t1, t_end;

  int thread_num_vec[3] = {1, 2, 4};

  gsl_rng ** generators = (gsl_rng **)malloc(omp_get_max_threads()*(sizeof(gsl_rng *)));
  for ( i = 0; i < omp_get_max_threads(); i++) generators[i] = gsl_rng_alloc(gsl_rng_taus);

  count = 0;
  thread_num = thread_num_vec[j];

  t0 = omp_get_wtime();
  #pragma omp parallel for reduction(+:count)
      for(i=0; i<trials; i++)
      {
        count += casino_game(generators[omp_get_thread_num()]);
      }
  t1 = omp_get_wtime();

  t_end = t1 - t0;


  mean_count = ((double) count )/( (double) trials );
  printf("number of threads: %d, end mean count: %f, end wall time: %f seconds \n", omp_get_max_threads(), mean_count, t_end);


  for ( i = 0; i < omp_get_max_threads(); i++) gsl_rng_free(generators[i]);
  free(generators);

  return 0;
}
