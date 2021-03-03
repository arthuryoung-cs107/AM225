#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <string.h>

#include "omp.h"
#include "Cash_Karp.hh"
#include "HW2_aux.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

int main()
{
  int i;

  double * lambda_vec = dvector(0, 11);
  int * tag_vec = ivector(0, 11);
  double evals, lambda_local;
  char prefix[100];
  memset(prefix, 0, 99);
  snprintf(prefix, 100, "./dat_dir/prob2_Bruss_eval_results");

  lambda_vec[0] = 1e-3;
  lambda_vec[1] = 1e-4;
  lambda_vec[2] = 1e-5;
  lambda_vec[3] = 1e-6;
  lambda_vec[4] = 1e-7;
  lambda_vec[5] = 1e-8;
  lambda_vec[6] = 1e-9;
  lambda_vec[7] = 1e-10;
  lambda_vec[8] = 1e-11;
  lambda_vec[9] = 1e-12;
  lambda_vec[10] = 1e-13;
  lambda_vec[11] = 1e-15;

  tag_vec[0] = -3;
  tag_vec[1] = -4;
  tag_vec[2] = -5;
  tag_vec[3] = -6;
  tag_vec[4] = -7;
  tag_vec[5] = -8;
  tag_vec[6] = -9;
  tag_vec[7] = -10;
  tag_vec[8] = -11;
  tag_vec[9] = -12;
  tag_vec[10] = -13;
  tag_vec[11] = -15;

  Brusselator * B_ref;
  FILE * eval_file = fopen("./dat_dir/prob2_Bruss_eval_results.aydat", "wb");
  for ( i = 0; i < 12; i++)
  {
    Brusselator * B1;
    B1 = new Brusselator(lambda_vec[i], lambda_vec[i], tag_vec[i]);
    B1->solve(0, 20);
    evals = (double) B1->evals;
    lambda_local = lambda_vec[i];
    fwrite(&(lambda_local), sizeof(double), 1, eval_file);
    fwrite(&(evals), sizeof(double), 1, eval_file);
    delete B1;
  }
  fclose(eval_file);
  aysml_gen(prefix, 12, 2);
  return 0;
}
