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
  // prob6_part_a();
  prob6_part_b();
  prob6_part_c();
  return 0;
}
