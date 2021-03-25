#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <string.h>

#include "omp.h"
#include "poisson_fft_AY.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

void prob5_part_a()
{
  int N = 112;
  int N_sqr = 21;
  char prefix[200];
  poisson_fft pf(N-1);

  memset(prefix, 0, 199);
  snprintf(prefix, 100, "./dat_dir/prob5_Poisson_fullgrid");

  pf.init_mms();

  pf.solve();
  pf.output_solution(prefix);
  printf("%d %g\n", N-1, pf.l2_error_mms() );
}

int main()
{
  prob5_part_a();
  return 0;
}
