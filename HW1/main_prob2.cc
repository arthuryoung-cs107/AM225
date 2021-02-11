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
  int i, j;
  double local_val;
  int generations = 500;
  char specfile[100];
  double * out_vec = dvector(0, 2);

  int n_vec[7] = {16, 32, 64, 128, 256, 512, 1024};

  memset(specfile, 0, 99);
  snprintf(specfile, 100, "./dat_dir/prob2_specs_threads%d.dat", omp_get_max_threads());
  FILE * prob2_specs_file = fopen(specfile, "w");
  fprintf(prob2_specs_file, "7 3");
  fclose(prob2_specs_file);

  memset(specfile, 0, 99);
  snprintf(specfile, 100, "./dat_dir/prob2_results_threads%d.dat", omp_get_max_threads());
  FILE * prob2_file = fopen(specfile, "wb");
  for ( i = 0; i < 7; i++)
  {
    cell_automation(n_vec[i], n_vec[i], generations, out_vec);
    fwrite(&(out_vec[0]), sizeof(double), 1, prob2_file);
    fwrite(&(out_vec[1]), sizeof(double), 1, prob2_file);
    fwrite(&(out_vec[2]), sizeof(double), 1, prob2_file);
  }
  fclose(prob2_file);

}
