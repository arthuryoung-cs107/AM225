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
  uint64_t i, j;
  uint64_t trials = 1e9;
  uint64_t count = 0;
  uint64_t jump = 100;
  int max_threads = 4;
  int thread_num;
  double mean_count, t0, t1, t_end;

  uint64_t start_seed_vec[4] = {888888888, 5555555555, 77777777777, 9999999999};
  int thread_num_vec[4] = {1, 2, 4};

  uint64_t * seed_vec = (uint64_t *) malloc(max_threads*sizeof(u_int64_t));

  for ( i = 0; i < max_threads; i++)
  {
    seed_vec[i] = lcg_fwd(start_seed_vec[i], jump);
  }

  for ( j = 0; j < 3; j++)
  {
    thread_num = thread_num_vec[j];
    t0 = omp_get_wtime();
    #pragma omp parallel for num_threads(thread_num)
        for(i=0; i<trials; i++)
        {
          uint64_t local_count = casino_game(seed_vec + omp_get_thread_num() );
    #pragma omp atomic
          count += local_count;
        }
    t1 = omp_get_wtime();

    t_end = t1 - t0;
    mean_count = ((double) count )/( (double) trials );
    printf("number of threads: %d, end mean count: %f, end wall time: %f seconds \n", thread_num, mean_count, t_end);
  }

  return 0;


}
