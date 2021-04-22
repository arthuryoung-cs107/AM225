#include "blas.h"
#include "omp.h"
#include "poisson_fft_AY.hh"
#include "Ritz_Galerk_sphere.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

void prob2_part_a()
{
  int N_full = 112;
  int N = N_full - 1;
  char prefix[200];
  char prefix2[200];

  auto f_main = [](double x, double y) { return std::exp(x-y); };

  poisson_fft pf(N, 1./((double) N+1));
  memset(prefix, 0, 199);
  memset(prefix2, 0, 199);
  snprintf(prefix, 100, "./dat_dir/prob2_square_fft");
  snprintf(prefix2, 100, "./dat_dir/prob2_square_fft_source");

  pf.init(f_main);
  pf.output_solution(prefix2, pf.f);
  pf.solve();
  pf.output_solution(prefix, pf.v);

  printf("entering Ritz Galerk init\n");
  Ritz_Galerk_sphere * FE1 = new Ritz_Galerk_sphere(N);

}

int main()
{

  prob2_part_a();

  // int m = 3; int n = 2;
  //
  // AYmat * A = new AYmat(m, n);
  // AYmat * B = new AYmat(n, 1);
  // AYmat * C = new AYmat(m, 1);
  //
  // A->init_123();
  // B->set(0, 0, 3.); B->set(1, 0, 4.);
  //
  // printf("A:\n");
  // A->print_mat();
  //
  // printf("B:\n");
  // B->print_mat();
  //
  // C->mult_set(A, B, 1., 0.);
  //
  // printf("C:\n");
  // C->print_mat();
  //
  // A->mult_put(B, C, 1., 0.);
  //
  // printf("C:\n");
  // C->print_mat();

  return 0;
}
