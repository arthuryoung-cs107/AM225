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
  uint64_t i, j, count, carry;
  uint64_t trials = 1e7;
  uint64_t jump = 10;
  uint64_t start_seed = 123456789;
  int thread_num = 2;
  double mean_count, t0, t1, t_end;

  uint64_t * seed_vec = (uint64_t *) malloc(thread_num*sizeof(uint64_t));

  for ( i = 0; i < thread_num; i++)
  {
    seed_vec[i] = lcg_fwd(start_seed+i, jump);
  }

  carry = lcg_fwd(start_seed, jump);
  count = 0;
  t0 = omp_get_wtime();
  #pragma omp parallel for num_threads(thread_num)
      for(i=0; i<trials; i++)
      {
        uint64_t local_count = casino_game(seed_vec + omp_get_thread_num() );
        count += local_count;
      }
  t1 = omp_get_wtime();

  t_end = t1 - t0;
  mean_count = ((double) count )/( (double) trials );
  printf("number of threads: %d, end mean count: %f, end wall time: %f seconds \n", thread_num, mean_count, t_end);

  return 0;


}
