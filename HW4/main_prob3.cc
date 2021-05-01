#include <cstring>

#include "blas.h"
#include "omp.h"
#include "cubic_1d_fe.hh"
#include "cubic_1d_alt_C2.hh"

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
  char prefix3[200];

  // N = 10;

  for ( N = 10; N <= 1000; N+=10)
  {
    cubic_1d_alt_C2 * FE1 = new cubic_1d_alt_C2(N);
    FE1->g = exp(-1.0)*(5.0)*M_PI;
    FE1->assemble_b();
    FE1->solve();
    memset(prefix, 0, 199);
    snprintf(prefix, 100, "./dat_dir/prob3_altcube_N%d_Ntest%d", N, N_test);
    FE1->write_out(prefix, N_test);
    delete FE1;
  }

  // //debugging zone
  N = 10;
  cubic_1d_alt_C2 * FE1 = new cubic_1d_alt_C2(N);
  FE1->assemble_b();
  double * omega = FE1->omega_1;
  double ** phi_check_all = dmatrix( 0, N_test-1, 0, 2 );
  double ** phi_check = dmatrix( 0, N_test-1, 0, FE1->n_full );
  double ** grad_phi_check = dmatrix( 0, N_test-1, 0, FE1->n_full );
  double del = (omega[1] - omega[0])/(N_test-1);
  double acc1, acc2;

  for ( i = 0; i < N_test; i++)
  {
    phi_check[i][0] = omega[0]+del*((double) i);
    grad_phi_check[i][0] = omega[0]+del*((double) i);
    phi_check_all[i][0] = omega[0]+del*((double) i);
  }

  memset(prefix, 0, 199);
  snprintf(prefix, 100, "./dat_dir/phi_check");
  memset(prefix2, 0, 199);
  snprintf(prefix2, 100, "./dat_dir/grad_phi_check");
  memset(prefix3, 0, 199);
  snprintf(prefix3, 100, "./dat_dir/phi_check_all");

  for ( i = 0; i < N_test; i++)
  {
    acc1 = 0.0;
    acc2 = 0.0;
    for ( j = 0; j < FE1->n_full; j++)
    {

      phi_check[i][j + 1] = FE1->phi_C2(phi_check[i][0], j) + ((double ) j);
      grad_phi_check[i][j + 1] = FE1->grad_phi_C2(phi_check[i][0], j) + ((double ) j);
      acc1 += FE1->phi_C2(phi_check[i][0], j);
      acc2 += FE1->grad_phi_C2(phi_check[i][0], j);
    }
    phi_check_all[i][1] = acc1;
    phi_check_all[i][2] = acc2;
  }
  fprintf_matrix(phi_check, N_test, FE1->n_full+1, prefix);
  fprintf_matrix(grad_phi_check, N_test, FE1->n_full+1, prefix2);
  fprintf_matrix(phi_check_all, N_test, 3, prefix3);

  // debugging zone

}

int main()
{
  prob3_part_a();

  return 0;
}
