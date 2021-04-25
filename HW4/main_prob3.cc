#include "blas.h"
#include "omp.h"
#include "cubic_1d_fe.hh"
#include "cubic_1d_alt.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

void prob3_part_a()
{
  int i,j;
  int N = 3;

  cubic_1d_alt * FE1 = new cubic_1d_alt(N);
}

int main()
{
  prob3_part_a();

  return 0;
}
