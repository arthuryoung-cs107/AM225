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
int ** Pierce_mat_gen(char pierce_filename[] , int m, int n )
{
  // note: we ignore edge and corner values, since they're outside the domain of the actual problem. Asusme
  int row_count, col_count, i, j;
  int ** pierce_mat = imatrix(1, m, 1, n);
  char extract[1000];
  zeromint_init(pierce_mat, 1, m, 1, n);

  FILE * pierce_file = fopen(pierce_filename, "r");

  while (fgets(extract, sizeof (extract), pierce_file))
  {
    row_count++;
    col_count = 0;
    for (i = 0; i < 1000; i++)
    {
      if (extract[i] == '0')
      {
        col_count++;
        pierce_mat[row_count][col_count] = 0;
        // printf("0 ");
      }
      else
      {
          if (extract[i] == '1')
          {
            col_count++;
            pierce_mat[row_count][col_count] = 1;
            // printf("1 ");
          }
      }
    }
    // printf("\n");
  }
  fclose(pierce_file);
  return pierce_mat;
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
