#include <cstring>

#include "blas.h"
#include "omp.h"
#include "cubic_1d_fe.hh"
#include "cubic_1d_alt_C2.hh"
#include "cubic_1d_alt_C1.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

void prob3_part_b_C1()
{
  int i,j, N;
  int N_test = 1000;
  double t0, t1, t_end;

  char prefix[200];
  char prefix2[200];
  char prefix3[200];
  char prefix4[200];

  double ** mat_out = dmatrix(0, 198, 0, 2);
  memset(prefix4, 0, 199);
  snprintf(prefix4, 100, "./dat_dir/prob3_altcube_C1_times");

  int count = 0;
  for ( N = 10; N <= 1000; N+=5)
  {
    cubic_1d_alt_C1 * FE1 = new cubic_1d_alt_C1(N);
    FE1->g = exp(-1.0)*(5.0)*M_PI;
    FE1->assemble_b();

    t0 = omp_get_wtime();
    FE1->solve();
    t1 = omp_get_wtime();

    memset(prefix, 0, 199);
    snprintf(prefix, 100, "./dat_dir/prob3_altcube_C1_N%d_Ntest%d", N, N_test);
    FE1->write_out(prefix, N_test);

    memset(prefix, 0, 199);
    snprintf(prefix, 100, "./dat_dir/prob3_altcube_C1_N%d_xsol", N);
    fprintf_matrix(&(FE1->x), 1, FE1->n, prefix4);


    t_end = t1 - t0;
    mat_out[count][0] = (double) N;
    mat_out[count][1] = FE1->h;
    mat_out[count][2] = t_end;

    delete FE1;
    count++;
  }
  fprintf_matrix(mat_out, count, 3, prefix4);


  //debugging zone
  N = 3;
  cubic_1d_alt_C1 * FE1 = new cubic_1d_alt_C1(N);
  FE1->assemble_b();
  double * omega = FE1->omega;
  double ** phi_check_all = dmatrix( 0, N_test-1, 0, 2 );
  double ** phi_check = dmatrix( 0, N_test-1, 0, FE1->n );
  double ** grad_phi_check = dmatrix( 0, N_test-1, 0, FE1->n );
  double del = (omega[1] - omega[0])/(N_test-1);
  double acc1, acc2;

  for ( i = 0; i < N_test; i++)
  {
    phi_check[i][0] = omega[0]+del*((double) i);
    grad_phi_check[i][0] = omega[0]+del*((double) i);
    phi_check_all[i][0] = omega[0]+del*((double) i);
  }

  memset(prefix, 0, 199);
  snprintf(prefix, 100, "./dat_dir/phi_check_C1");
  memset(prefix2, 0, 199);
  snprintf(prefix2, 100, "./dat_dir/grad_phi_check_C1");
  memset(prefix3, 0, 199);
  snprintf(prefix3, 100, "./dat_dir/phi_check_all_C1");

  for ( i = 0; i < N_test; i++)
  {
    acc1 = 0.0;
    acc2 = 0.0;
    for ( j = 0; j < FE1->n; j++)
    {

      phi_check[i][j + 1] = FE1->phi_C1(phi_check[i][0], j) + ((double ) j);
      grad_phi_check[i][j + 1] = FE1->grad_phi_C1(phi_check[i][0], j) + ((double ) j);
      acc1 += FE1->phi_C1(phi_check[i][0], j);
      acc2 += FE1->grad_phi_C1(phi_check[i][0], j);
    }
    phi_check_all[i][1] = acc1;
    phi_check_all[i][2] = acc2;
  }
  fprintf_matrix(phi_check, N_test, FE1->n+1, prefix);
  fprintf_matrix(grad_phi_check, N_test, FE1->n+1, prefix2);
  fprintf_matrix(phi_check_all, N_test, 3, prefix3);

  // debugging zone

}

void prob3_part_b_C2()
{
  int i,j, N;
  int N_test = 1000;

  double t0, t1, t_end;

  char prefix[200];
  char prefix2[200];
  char prefix3[200];
  char prefix4[200];

  double ** mat_out = dmatrix(0, 198, 0, 2);
  memset(prefix4, 0, 199);
  snprintf(prefix4, 100, "./dat_dir/prob3_altcube_C2_times");

  // N = 10;
  int count = 0;
  for ( N = 10; N <= 1000; N+=5)
  {
    cubic_1d_alt_C2 * FE1 = new cubic_1d_alt_C2(N);
    FE1->g = exp(-1.0)*(5.0)*M_PI;
    FE1->assemble_b();

    t0 = omp_get_wtime();
    FE1->solve();
    t1 = omp_get_wtime();

    memset(prefix, 0, 199);
    snprintf(prefix, 100, "./dat_dir/prob3_altcube_C2_N%d_Ntest%d", N, N_test);
    FE1->write_out(prefix, N_test);

    memset(prefix, 0, 199);
    snprintf(prefix, 100, "./dat_dir/prob3_altcube_C2_N%d_xsol", N);
    fprintf_matrix(&(FE1->x), 1, FE1->n, prefix4);

    t_end = t1 - t0;
    mat_out[count][0] = (double) N;
    mat_out[count][1] = FE1->h;
    mat_out[count][2] = t_end;

    delete FE1;
    count++;
  }

  fprintf_matrix(mat_out, count, 3, prefix4);

  // //debugging zone
  N = 3;
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
  snprintf(prefix, 100, "./dat_dir/phi_check_C2");
  memset(prefix2, 0, 199);
  snprintf(prefix2, 100, "./dat_dir/grad_phi_check_C2");
  memset(prefix3, 0, 199);
  snprintf(prefix3, 100, "./dat_dir/phi_check_all_C2");

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

void prob3_part_b_C0()
{
  int i_max = 30;
  double mean_count, t0, t1, t_end;
  double ** output_mat = dmatrix(0, i_max, 0, 3);

  char prefix[200];
  memset(prefix, 0, 199);
  snprintf(prefix, 100, "./dat_dir/prob3_cube_C0_error");

  for(int i=0;i<=i_max;i++)
  {
      // Create the finite-element problem and set the Neumann boundary
      // condition to match the manufactured solution
      int j=int(10*pow(100,(1/30.)*i)+0.5);
      cubic_1d_fe cf(j);
      cf.g=exp(-1)*5*M_PI;

      // Initialize the source term for the manufactured solution, solve, and
      // the print the L2 error
      cf.init_mms();

      t0 = omp_get_wtime();
      cf.solve();
      t1 = omp_get_wtime();
      t_end = t1-t0;

      output_mat[i][0] = (double) j;
      output_mat[i][1] = cf.h;
      output_mat[i][2] = t_end;
      output_mat[i][3] = cf.l2_norm_mms();

      printf("%d %g %g\n",j,cf.h,cf.l2_norm_mms());
  }
  fprintf_matrix(output_mat, i_max, 4, prefix);

}

int main()
{
  prob3_part_b_C2();
  prob3_part_b_C1();
  prob3_part_b_C0();

  return 0;
}
