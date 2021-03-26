#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <string.h>

#include "omp.h"
#include "poisson_fft_AY.hh"
#include "square_specs.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}
void prob5_part_a()
{
  int i, j;
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

  square_specs * S1 = new square_specs();

}

int main()
{
  int i, j;

  const char * check1 = ".";
  const char * check2 = "1";
  const char * check3 = "U";
  printf("%s : %d \n", check1, (int) *check1);
  printf("%s : %d \n", check2, (int) *check2);
  printf("%s : %d \n", check3, (int) *check3);

  char ascii_name[22][2] = {".", "A", "B", "C", "D", "E", "F", "G", "H", "I",
  "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U"};
  int ascii_num[22];
  for ( i = 0; i < 22; i++) ascii_num[i]= (int) *(ascii_name[i]);

  printf("%s : %d \n", ascii_name[0], (int) ascii_name[0][0]);
  printf("%s : %d, %c \n", ascii_name[1],  ascii_num[1], ascii_num[1]);

  prob5_part_a();
  return 0;
}
