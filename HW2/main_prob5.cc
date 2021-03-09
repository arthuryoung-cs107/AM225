#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <string.h>

#include "omp.h"
#include "HW2_aux.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

int main()
{
  // prob5_part_a();
  prob5_part_b();
  return 0;
}
