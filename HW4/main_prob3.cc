#include <cstring>

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
  int i,j, N;
  int N_test = 1000;

  char prefix[200];
  char prefix2[200];

  // N = 100;

  // for ( N = 10; N <= 1000; N+=10)
  // {
  //   cubic_1d_alt * FE1 = new cubic_1d_alt(N);
  //   FE1->g = exp(-1.0)*(5.0)*M_PI;
  //   FE1->assemble_b();
  //   FE1->solve();
  //
  //   memset(prefix, 0, 199);
  //   snprintf(prefix, 100, "./dat_dir/prob3_altcube_N%d_Ntest%d", N, N_test);
  //
  //   FE1->write_out(prefix, N_test);
  //   delete FE1;
  // }

  //debugging zone
  N = 4;
  double ** phi_check = dmatrix( 0, N_test-1, 0, N );
  double ** grad_phi_check = dmatrix( 0, N_test-1, 0, N );
  double del = (1.0)/(N_test-1);

  for ( i = 0; i < N_test; i++)
  {
    phi_check[i][0] = 1.0+del*((double) i);
    grad_phi_check[i][0] = 1.0+del*((double) i);
  }

  cubic_1d_alt * FE1 = new cubic_1d_alt(N);
  memset(prefix, 0, 199);
  snprintf(prefix, 100, "./dat_dir/phi_check");
  memset(prefix2, 0, 199);
  snprintf(prefix2, 100, "./dat_dir/grad_phi_check");

  for ( i = 0; i < N_test; i++)
  {
    for ( j = 0; j < N; j++)
    {
      phi_check[i][j + 1] = FE1->phi_C1(phi_check[i][0], j) + ((double ) j);
      grad_phi_check[i][j + 1] = FE1->grad_phi_C1(phi_check[i][0], j) + ((double ) j);
    }
  }
  fprintf_matrix(phi_check, N_test, N+1, prefix);
  fprintf_matrix(grad_phi_check, N_test, N+1, prefix2);

  //debugging zone

}

int main()
{
  prob3_part_a();

  return 0;
}
