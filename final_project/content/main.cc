#include "blas.h"
#include "omp.h"
#include "ADM.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

void test1()
{
  int m = 4; int n = 3;
  AYmat * mat1 = new AYmat(m, n);
  AYmat * mat2 = new AYmat(m, n);
  AYmat * mat3 = new AYmat(m, n);
  AYmat * mat4 = new AYmat(m, n);

  mat1->init_123();
  mat2->init_123();
  mat3->init_123();

  mat1->add(mat2, mat3, -2., 0.);
  printf("\nmat3\n");
  mat3->print_mat();

  mat3->add(mat2, mat3, 2., 1.);
  printf("\nmat3\n");
  mat3->print_mat();

  mat1->copy_set(mat4);
  printf("\nmat4\n");
  mat4->print_mat();

}

void test2()
{
  double t0, t1, tend;
  char specfile_X_noised[200];
  memset(specfile_X_noised, 0, 199);
  snprintf(specfile_X_noised, 200, "./logo_dat_ldir/X_noised");

  char prefix1[200];
  memset(prefix1, 0, 199);
  snprintf(prefix1, 200, "./logo_dat_dir/test2");

  AYmat * X_noised_main = aysml_read(specfile_X_noised);

  ADM * adm1 = new ADM(X_noised_main);

  t0 = omp_get_wtime();
  adm1->solve(true);
  t1 = omp_get_wtime();

  adm1->write_out(prefix1);

  delete adm1;

  tend = t1-t0;
  printf("time elapsed: %f seconds\n", tend);
}

void test3()
{
#pragma omp parallel for
  for (int i = 500; i >= 0; i--)
  {

    double t0, t1, tend;
    char specfile_X_noised[200];
    memset(specfile_X_noised, 0, 199);
    snprintf(specfile_X_noised, 200, "./aydat_dir_small_sim3/tem%d", i );

    char prefix1[200];
    memset(prefix1, 0, 199);
    snprintf(prefix1, 200, "./aydat_dir_small_sim3/sim1_%d", i);

    AYmat * X_noised_main_T = aysml_read(specfile_X_noised);
    AYmat * X_noised_main = X_noised_main_T->transpose_gen();

    delete X_noised_main_T;

    ADM * adm1 = new ADM(X_noised_main);
    delete X_noised_main;

    t0 = omp_get_wtime();
    adm1->solve();
    t1 = omp_get_wtime();

    adm1->write_out(prefix1);

    delete adm1;

    tend = t1-t0;
    printf("field %d done. Time elapsed: %f seconds\n", i, tend);

  }

}

void test3_corrupt()
{
#pragma omp parallel for
  for (int i = 500; i >= 0; i--)
  {

    double t0, t1, tend;
    char specfile_X_noised[200];
    memset(specfile_X_noised, 0, 199);
    snprintf(specfile_X_noised, 200, "./aydat_dir_small_sim3/tem_corrupt%d", i );

    char prefix1[200];
    memset(prefix1, 0, 199);
    snprintf(prefix1, 200, "./aydat_dir_small_sim3/sim_corrupt_%d", i);

    AYmat * X_noised_main_T = aysml_read(specfile_X_noised);
    AYmat * X_noised_main = X_noised_main_T->transpose_gen();

    delete X_noised_main_T;

    ADM * adm1 = new ADM(X_noised_main);
    delete X_noised_main;

    t0 = omp_get_wtime();
    adm1->solve();
    t1 = omp_get_wtime();

    adm1->write_out(prefix1);

    delete adm1;

    tend = t1-t0;
    printf("field %d done. Time elapsed: %f seconds\n", i, tend);

  }

}

void test4()
{
#pragma omp parallel for
  for (int i = 500; i >= 0; i--)
  {

    double t0, t1, tend;
    char specfile_X_noised[200];
    memset(specfile_X_noised, 0, 199);
    snprintf(specfile_X_noised, 200, "./aydat_dir_small_sim4/tem%d", i );

    char prefix1[200];
    memset(prefix1, 0, 199);
    snprintf(prefix1, 200, "./aydat_dir_small_sim4/sim1_%d", i);

    AYmat * X_noised_main_T = aysml_read(specfile_X_noised);
    AYmat * X_noised_main = X_noised_main_T->transpose_gen();

    delete X_noised_main_T;

    ADM * adm1 = new ADM(X_noised_main);
    delete X_noised_main;

    t0 = omp_get_wtime();
    adm1->solve();
    t1 = omp_get_wtime();

    adm1->write_out(prefix1);

    delete adm1;

    tend = t1-t0;
    printf("field %d done. Time elapsed: %f seconds\n", i, tend);

  }

}

void test4_corrupt()
{
#pragma omp parallel for
  for (int i = 500; i >= 0; i--)
  {

    double t0, t1, tend;
    char specfile_X_noised[200];
    memset(specfile_X_noised, 0, 199);
    snprintf(specfile_X_noised, 200, "./aydat_dir_small_sim4/tem_corrupt%d", i );

    char prefix1[200];
    memset(prefix1, 0, 199);
    snprintf(prefix1, 200, "./aydat_dir_small_sim4/sim_corrupt_%d", i);

    AYmat * X_noised_main_T = aysml_read(specfile_X_noised);
    AYmat * X_noised_main = X_noised_main_T->transpose_gen();

    delete X_noised_main_T;

    ADM * adm1 = new ADM(X_noised_main);
    delete X_noised_main;

    t0 = omp_get_wtime();
    adm1->solve();
    t1 = omp_get_wtime();

    adm1->write_out(prefix1);

    delete adm1;

    tend = t1-t0;
    printf("field %d done. Time elapsed: %f seconds\n", i, tend);

  }

}

void test5()
{
#pragma omp parallel for
  for (int i = 500; i >= 0; i--)
  {

    double t0, t1, tend;
    char specfile_X_noised[200];
    memset(specfile_X_noised, 0, 199);
    snprintf(specfile_X_noised, 200, "./aydat_dir_small_sim5/tem%d", i );

    char prefix1[200];
    memset(prefix1, 0, 199);
    snprintf(prefix1, 200, "./aydat_dir_small_sim5/sim1_%d", i);

    AYmat * X_noised_main_T = aysml_read(specfile_X_noised);
    AYmat * X_noised_main = X_noised_main_T->transpose_gen();

    delete X_noised_main_T;

    ADM * adm1 = new ADM(X_noised_main);
    delete X_noised_main;

    t0 = omp_get_wtime();
    adm1->solve();
    t1 = omp_get_wtime();

    adm1->write_out(prefix1);

    delete adm1;

    tend = t1-t0;
    printf("field %d done. Time elapsed: %f seconds\n", i, tend);

  }

}

void test5_corrupt()
{
#pragma omp parallel for
  for (int i = 500; i >= 0; i--)
  {

    double t0, t1, tend;
    char specfile_X_noised[200];
    memset(specfile_X_noised, 0, 199);
    snprintf(specfile_X_noised, 200, "./aydat_dir_small_sim5/tem_corrupt%d", i );

    char prefix1[200];
    memset(prefix1, 0, 199);
    snprintf(prefix1, 200, "./aydat_dir_small_sim5/sim_corrupt_%d", i);

    AYmat * X_noised_main_T = aysml_read(specfile_X_noised);
    AYmat * X_noised_main = X_noised_main_T->transpose_gen();

    delete X_noised_main_T;

    ADM * adm1 = new ADM(X_noised_main);
    delete X_noised_main;

    t0 = omp_get_wtime();
    adm1->solve();
    t1 = omp_get_wtime();

    adm1->write_out(prefix1);

    delete adm1;

    tend = t1-t0;
    printf("field %d done. Time elapsed: %f seconds\n", i, tend);

  }

}

int main()
{
  // test1();
  // test2();
  // test3();
  // test4();
  // test5();

  // test3_corrupt();
  // test4_corrupt();
  // test5_corrupt();


  return 0;
}
