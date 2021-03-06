#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdint.h>
#include <string.h>

#include "blas.h"
#include "omp.h"
#include "poisson_fft_AY.hh"
#include "square_specs.hh"
#include "schur_perfect.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}
void prob5_part_a()
{
  int i, j;
  int N_full = 112;
  int N_sqr = 21;
  char prefix[200];
  char prefix2[200];

  int N = N_full-1;
  auto f_main = [](double x, double y) { return std::exp(x-y); };

  poisson_fft pf(N, 1./((double) N+1));
  memset(prefix, 0, 199);
  memset(prefix2, 0, 199);
  snprintf(prefix, 100, "./dat_dir/prob5_Poisson_fullgrid");
  snprintf(prefix2, 100, "./dat_dir/prob5_Poisson_fullgrid_source");

  pf.init(f_main);
  pf.output_solution(prefix2, pf.f);
  pf.solve();
  pf.output_solution(prefix, pf.v);
}
void prob5_part_b()
{
  int i, j, k, l ;
  char prefix[200];
  char prefix2[200];

  memset(prefix, 0, 199);
  memset(prefix2, 0, 199);
  snprintf(prefix, 100, "./dat_dir/prob5_Poisson_Schur");
  snprintf(prefix2, 100, "./dat_dir/prob5_Poisson_Schur_source");

  square_specs * S1 = new square_specs();
  // S1->print_adjacents();

  auto f = [](double x, double y) { return std::exp(x-y); };

  schur_perfect * P1 = new schur_perfect(S1);
  P1->init(f);
  P1->output_solution(prefix2, P1->f_full);
  P1->solve_S();
  P1->output_solution(prefix, P1->v_sol);

}
int main()
{

  prob5_part_a();
  prob5_part_b();

  return 0;
}
