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
  Brusselator B1;
  B1.solve(0, 20);
  return 0;
}
