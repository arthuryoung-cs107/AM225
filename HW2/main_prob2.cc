#include <cstdio>
#include <cstdlib>
#include <stdint.h>

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

  double * lambda_vec = dvector(0, 10);
  int * tag_vec = ivector(0, 10);
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

  Brusselator * B_ref;
  // B_ref(lambda_vec[i], lambda_vec[i], tag_vec[i]);


  for ( i = 0; i < 11; i++)
  {
    Brusselator * B1;
    B1 = new Brusselator(lambda_vec[i], lambda_vec[i], tag_vec[i]);
    B1->solve(0, 20);
    delete B1;
  }

  return 0;
}
