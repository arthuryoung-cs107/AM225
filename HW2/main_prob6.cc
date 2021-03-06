#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <string.h>

#include "omp.h"
#include "Geng.hh"
#include "HW2_aux.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

int main()
{
  Brusselator_Geng * Brus_Geng;
  Brus_Geng = new Brusselator_Geng();

  printf("survived\n");

  return 0;
}
