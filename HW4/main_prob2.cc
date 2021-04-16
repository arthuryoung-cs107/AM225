#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdint.h>
#include <string.h>

#include "blas.h"
#include "omp.h"
#include "poisson_fft_AY.hh"

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

  Ritz_Galerk * FE1 = new Ritz_Galerk();
  // sphere_Ritz_Galerk * FE1 = new sphere_Ritz_Galerk(N);

}

int main()
{

  prob2_part_a();

  return 0;
}
