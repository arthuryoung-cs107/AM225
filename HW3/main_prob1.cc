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
  int i, j, n, trials;
  trials = 4;
  for ( n = 16; n <= 4096; n *= 2)
  {
    double net_time = 0;
    #pragma omp parallel for reduction(+:net_time)

      for ( i = 0; i < trials; i++)
      {
        AYmat * Mat1 = new AYmat(n, n);
        AYmat * Mat2 = new AYmat(n, n);
        AYmat * Mat3 = new AYmat(n, n);

        Mat1->init_randuni();
        Mat2->init_randuni();
        Mat3->init_randuni();
        double t0 = omp_get_wtime();
        AYmat_mul(Mat1, Mat2, Mat3);
        double t1 = omp_get_wtime();
        net_time += (t1 - t0);

        delete Mat1;
        delete Mat2;
        delete Mat3;
      }
      double t_avg = net_time/((double) trials );
    printf("%d %e\n", n, t_avg );

  }

}

void prob1_part_b()
{
  int i, j, n, trials;
  trials = 4;

    AYmat * Mat1 = new AYmat(16, 16);
    AYmat * Mat2 = new AYmat(16, 16);
    AYmat * Mat3 = new AYmat(16, 16);

    Mat1->init_randuni();
    Mat2->init_randuni();

    Mat1->print_mat();
    printf("\n");
    Mat2->print_mat();
    printf("\n");

    AYmat_mul_Strass(Mat1, Mat2, Mat3);

    Mat3->print_mat();
    
  // for ( n = 16; n < 4096; n *= 2)
  // {
  //   AYmat * Mat1 = new AYmat(n, n);
  //   AYmat * Mat2 = new AYmat(n, n);
  //   AYmat * Mat3 = new AYmat(n, n);
  //
  //   Mat1->init_randuni();
  //   Mat2->init_randuni();
  //   Mat3->init_randuni();
  //
  //   Mat1
  //
  //   AYmat_mul_Strass(Mat1, Mat2, Mat3);
  //
  //   delete Mat1;
  //   delete Mat2;
  //   delete Mat3;
  // }
}
int main()
{
  // prob1_part_a();
  prob1_part_b();

  return 0;
}
