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
  uint64_t i, j, count, carry, thread_num;
  uint64_t trials = 1e9;
  uint64_t jump = 10;
  uint64_t start_seed = 123456789;
  int max_threads = 4;
  double mean_count, t0, t1, t_end;

  uint64_t * seed_vec = (uint64_t *) malloc(max_threads*sizeof(uint64_t));
  int thread_num_vec[3] = {1, 2, 4};
  uint64_t count_vec[4] = {0, 0, 0, 0};

  for ( i = 0; i < max_threads; i++)
  {
    seed_vec[i] = lcg_fwd(start_seed+i, jump);
  }

  carry = lcg_fwd(start_seed, jump);

  for (size_t j = 0; j < 3; j++)
  {
    count_vec[0] = 0;
    count_vec[1] = 0;
    count_vec[2] = 0;
    count_vec[3] = 0;

    thread_num = thread_num_vec[j];

    t0 = omp_get_wtime();
    #pragma omp parallel for num_threads(thread_num)
        for(i=0; i<trials; i++)
        {
          uint64_t local_count = casino_game(seed_vec + omp_get_thread_num() );
          count_vec[omp_get_thread_num()] += local_count;
        }
    t1 = omp_get_wtime();

    t_end = t1 - t0;

    count = count_vec[0] + count_vec[1] + count_vec[2] + count_vec[3];

    mean_count = ((double) count )/( (double) trials );
    printf("number of threads: %ld, end mean count: %f, end wall time: %f seconds \n", thread_num, mean_count, t_end);

  }

  return 0;


}
