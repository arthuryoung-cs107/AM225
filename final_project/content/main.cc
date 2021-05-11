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
  AYmat * mat2 = new AYmat(n, m);
  AYmat * mat3 = new AYmat(n, 1);

  mat1->init_123();
  mat2->init_123();
  mat3->init_123();

  // mat1->print_mat();
  // mat2->print_mat();

  AYmat * mat4 = mat1->mult_gen(mat2, 1.0);
  mat4->print_mat();

  mat1->mult_set(mat2, mat4, 1.0, 0.0);
  mat4->print_mat();

  AYmat * mat5 = mat1->mult_gen(mat3, 1.0);
  mat5->print_mat();

  mat1->mult_set(mat3, mat5, 1.0, 0.0);
  mat5->print_mat();
}

void test2()
{
  int m = 4; int n = 3;

  AYmat * X_main = new AYmat(m, n);
  X_main->init_randuni();
  X_main->print_mat();

  ADM * adm1 = new ADM(X_main);

  delete adm1;

}

int main()
{
  // test1();
  test2();

  return 0;
}
