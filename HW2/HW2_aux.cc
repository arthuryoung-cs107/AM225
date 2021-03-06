#include "HW2_aux.hh"
#include "omp.h"
#include "schemes.hh"

extern "C" {
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

void prob5_part_b()
{
  double * J_main = dvector(0, 1);
  double * K_main = dvector(0, 1);

  J_main[0] = 0.5;
  J_main[1] = 1.0;
  K_main[0] = 0.5;
  K_main[1] = -0.2;

  Kuramoto2D * K5, *K6;
  K5 = new Kuramoto2D(1e-6, 1e-6, J_main[0], K_main[0], 5);
  K5->solve(0, 200, 401);
  delete K5;

  K6 = new Kuramoto2D(1e-6, 1e-6, J_main[1], K_main[1], 6);
  K6->solve(0, 200, 401);
  delete K6;

  Kuramoto2D_SUPER * K56;
  K56 = new Kuramoto2D_SUPER(1e-6, 1e-6, J_main, K_main, 5, 6);
  K56->solve(0, 200, 401);
  delete K56;


}

void prob5_part_a()
{
  double J_vec[3] = {0.5, 0.3, 1.0};
  double K_vec[3] = {0.5, -0.2, -0.2};

  # pragma omp parrallel for
    for (int i = 0; i < 3; i++)
    {
      Kuramoto2D * K1;
      K1 = new Kuramoto2D(1e-6, 1e-6, J_vec[i], K_vec[i], i+1);
      K1->solve(0, 200, 401);
      delete K1;
    }
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

  free_ivector(step_vector, 1, 11);
}

void prob6_part_b()
{
  double * params_main = dvector(0, 5);

  params_main[0] = 1.0 ;
  params_main[1] = 0.25 ;
  params_main[2] = 1.0 ;
  params_main[3] = 1.25 ;
  params_main[4] = 1.0 ;
  params_main[5] = 0.75 ;

  Galaxy_Geng * G1;
  G1 = new Galaxy_Geng(params_main, 1);
  G1->solve(2000, 1);
  delete G1;

}

void prob6_part_c()
{
  double * params_main = dvector(0, 5);

  params_main[0] = 1.0 ;
  params_main[1] = 0.25 ;
  params_main[2] = 1.0 ;
  params_main[3] = 1.25 ;
  params_main[4] = 1.0 ;
  params_main[5] = 0.75 ;

  Galaxy_Geng * G2;
  G2 = new Galaxy_Geng(params_main, 2);
  G2->solve(1e5, 2);

  delete G2;
}
