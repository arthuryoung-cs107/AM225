#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdint.h>
#include <string.h>

#include "blas.h"
#include "omp.h"
#include "AYmat.hh"
#include "AYmat_ops.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}
void prob1_part_a()
{
  int i, j, n;

  for ( n = 16; n < 1024; n *= 2)
  {
    AYmat * Mat1 = new AYmat(n, n);
    AYmat * Mat2 = new AYmat(n, n);
    AYmat * Mat3 = new AYmat(n, n);

    Mat1->init_randuni();
    Mat2->init_randuni();
    Mat3->init_randuni();

    AYmat_mul(Mat1, Mat2, Mat3);

    delete Mat1;
    delete Mat2;
    delete Mat3;
  }

}

int main()
{
  prob1_part_a();

  return 0;
}
