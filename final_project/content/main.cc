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
  snprintf(specfile_X_noised, 200, "./dat_dir/X_noised");

  char prefix1[200];
  memset(prefix1, 0, 199);
  snprintf(prefix1, 200, "./dat_dir/test2");

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

int main()
{
  // test1();
  // test2();


  return 0;
}
