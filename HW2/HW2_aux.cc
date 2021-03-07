#include "HW2_aux.hh"
#include "omp.h"
#include "schemes.hh"


extern "C" {
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

void prob6_part_a()
{
  int i;
  // int * step_vector = ivector(1, 13);
  int * step_vector = ivector(1, 11);

  step_vector[1] = 1e2;
  step_vector[2] = 5e2;
  step_vector[3] = 1e3;
  step_vector[4] = 5e3;
  step_vector[5] = 1e4;
  step_vector[6] = 5e4;
  step_vector[7] = 1e5;
  step_vector[8] = 5e5;
  step_vector[9] = 1e6;
  step_vector[10] = 5e6;
  step_vector[11] = 1e7;
  // step_vector[12] = 5e7;
  // step_vector[13] = 1e8;

  # pragma omp parrallel for
    // for ( i = 1; i <= 13; i++)
    for ( i = 1; i <= 11; i++)
    {

      int tag1 = 1;
      int tag2;
      if (i%2 == 0)
      {
        tag1 = 5;
      }
      tag2 = (i+3)/((int) 2);

      Brusselator_Geng * Brus_Geng;
      Brus_Geng = new Brusselator_Geng(step_vector[i], tag1, tag2);
      Brus_Geng->solve(0.0, 20.0);
      delete Brus_Geng;
    }

}
