#include <cstring>

#include "blas.h"
#include "omp.h"
extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

void prob3_part_a()
{
  int N;
  int N_test = 40;

  char prefix_x[200];
  char prefix_y[200];
  char prefix_v[200];
  char prefix_w[200];
  char prefix_u[200];

  auto f_main = [](double v, double w) { return ((std::exp(-v))*( 3.0 + (v - 4.0)*v + w*w)); };

  for ( N = 10; N <= 100; N+=5)
  {
      memset(prefix_x, 0, 199);
      snprintf(prefix_x, 100, "./dat_dir/prob2_Ritz_Galerk_x_N%d_Ntest%d", N, N_test);
      memset(prefix_y, 0, 199);
      snprintf(prefix_y, 100, "./dat_dir/prob2_Ritz_Galerk_y_N%d_Ntest%d", N, N_test);
      memset(prefix_v, 0, 199);
      snprintf(prefix_v, 100, "./dat_dir/prob2_Ritz_Galerk_v_N%d_Ntest%d", N, N_test);
      memset(prefix_w, 0, 199);
      snprintf(prefix_w, 100, "./dat_dir/prob2_Ritz_Galerk_w_N%d_Ntest%d", N, N_test);
      memset(prefix_u, 0, 199);
      snprintf(prefix_u, 100, "./dat_dir/prob2_Ritz_Galerk_u_N%d_Ntest%d", N, N_test);

      Ritz_Galerk_sphere * FE1 = new Ritz_Galerk_sphere(N);
      FE1->assemble_b(f_main);
      FE1->solve();

      FE1->write_out(prefix_x, prefix_y, prefix_v, prefix_w, prefix_u, N_test);
      delete FE1;
  }

}

int main()
{

  prob3_part_a();

  return 0;
}
