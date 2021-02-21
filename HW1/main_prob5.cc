#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <string.h>

#include "HW1_aux.hh"
#include "omp.h"

extern "C" {
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

int main()
{
int i;
double runtime, thread_count_cast, N_cast;
int * N_vec = ivector(0, 5);
double ** T1 = dmatrix(0, 1, 0, 1);
double ** Tid = dmatrix(0, 1, 0, 1);
char prefix[100];
char specfile[200];

N_vec[0] = 64;
N_vec[1] = 128;
N_vec[2] = 256;
N_vec[3] = 512;
N_vec[4] = 1024;
N_vec[5] = 2048;

memset(prefix, 0, 99);
memset(specfile, 0, 199);
snprintf(prefix, 100, "./dat_dir/prob5_threads%d_runtime", omp_get_max_threads());
snprintf(specfile, 200, "%s.aydat", prefix);
FILE * prob5_data_file = fopen(specfile, "wb");

for ( i = 0; i < 6; i++)
{
  runtime = solve_grid(N_vec[i]);
  N_cast = (double) N_vec[i];
  thread_count_cast = (double) omp_get_max_threads();
  fwrite(&(thread_count_cast), sizeof(double), 1, prob5_data_file);
  fwrite(&(N_cast), sizeof(double), 1, prob5_data_file);
  fwrite(&(runtime), sizeof(double), 1, prob5_data_file);
}
fclose(prob5_data_file);

aysml_gen(prefix, 6, 3);

T1[0][0] = 1;
T1[0][1] = 0.5;
T1[1][0] = 0;
T1[1][1] = 1;

Tid[0][0] = 1;
Tid[0][1] = 0;
Tid[1][0] = 0;
Tid[1][1] = 1;

memset(prefix, 0, 99);
snprintf(prefix, 100, "./dat_dir/prob5_threads%d_s_vs_t_T1", omp_get_max_threads());
solve_grid_integrate(400, T1, prefix);

memset(prefix, 0, 99);
snprintf(prefix, 100, "./dat_dir/prob5_threads%d_s_vs_t_Tid", omp_get_max_threads());
solve_grid_integrate(400, Tid, prefix);

}
