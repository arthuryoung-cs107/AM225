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
  int N = 20;
  char prefix[200];
  char prefix2[200];

  auto f_main = [](double v, double w) { return ((std::exp(-v))*( 3.0 + (v - 4.0)*v + w*w)); };

  memset(prefix, 0, 199);
  snprintf(prefix, 100, "./dat_dir/prob2_Ritz_Galerk_xyvw");
  memset(prefix2, 0, 199);
  snprintf(prefix2, 100, "./dat_dir/prob2_Ritz_Galerk_u");

  Ritz_Galerk_sphere * FE1 = new Ritz_Galerk_sphere(N);
  FE1->assemble_b(f_main);
  FE1->solve(true);

  FE1->write_out(prefix, prefix2, 40);

}

int main()
{

  prob2_part_a();

  return 0;
}
