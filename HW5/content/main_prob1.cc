#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "blas.h"
#include "omp.h"
#include "fluid_2d.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

void prob1_part_a()
{
  int N;

  char prefix[200];
  memset(prefix, 0, 199);
  snprintf(prefix, 100, "./dat_dir/prob1");




}

int main()
{

  prob1_part_a();

  return 0;
}
