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
  // int trials = 1e9;
  int trials = 5;
  int thread_num = 1;
  uint64_t start_seed = 888888888;
  uint64_t jump = 100;

  uint64_t * seed_vec = (uint64_t *) malloc(thread_num*sizeof(uint64_t));
  int * count_vec = ivector(0, trials-1);

  seed_vec[0] = lcg_fwd(start_seed, jump);
  for ( i = 1; i < thread_num; i++)
  {
    start_seed += 1;
    seed_vec[i] = lcg_fwd(start_seed, jump);
  }

  #pragma omp parallel for num_threads(thread_num)
      for(i=0; i<trials; i++)
      {
        count_vec[i] = casino_game(seed_vec + omp_get_thread_num() );
        printf("trial: %d, thread: %d, count: %d \n", i, omp_get_thread_num(), count_vec[i]);
      }


  return 0;
}
