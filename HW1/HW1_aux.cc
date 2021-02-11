#include "HW1_aux.hh"
#include "omp.h"
#include "gsl/gsl_rng.h"

extern "C" {
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

double random_uni(double low, double high, uint64_t * carry) // knuth rng
{
  double rand_uni = ((double) lcg_uni(carry))/(lcg_sze());
  rand_uni = rand_uni*(high - low) + ((high - low)/2 - 0.5);
  return rand_uni;
}

uint64_t casino_game(gsl_rng * T)
{
  double end_cond = 1;
  double cond = 0;
  uint64_t count = 0;
  while(cond < end_cond)
  {
    cond += gsl_rng_uniform(T);
    count++;
  }
  return count;
}

int sum_range(int ** cell, int row_low, int row_high, int col_low, int col_high)
{
  int sum = 0;

  for (int i = row_low; i <= row_high; i++)
  {
    for (int j = col_low; j <= col_high; j++)
    {
      sum += cell[i][j];
    }
  }
  return sum;
}

int update_rule(int status, int sum)
{
  int s_surv = 5;
  int s_birth = 3;
  int return_val;

  if (status == 1)
  {
    if ( (sum >= 1) && (sum <= s_surv))
    {
      return_val = 1;
    }
    else
    {
      return_val = 0;
    }
  }
  else
  {
    if ( sum == s_birth )
    {
      return_val = 1;
    }
    else
    {
      return_val = 0;
    }
  }
  return return_val;
}

void update_row(int ** cell_old, int i, int m, int n, int ** cell)
{
  int j;
  int local_sum;
  if (i == 0)
  {
    cell[i][0] = update_rule(cell_old[i][0], sum_range(cell_old, i, i+1, 0, 1) - cell_old[i][0]);
    for ( j = 1; j < (n-1); j++)
    {
      cell[i][j] = update_rule(cell_old[i][j], sum_range(cell_old, i, i+1, j-1, j+1) - cell_old[i][j]);
    }
    cell[i][n-1] = update_rule(cell_old[i][n-1], sum_range(cell_old, i, i+1, n-2, n-1) - cell_old[i][n-1]);
  }
  else
  {
    if (i == (m-1))
    {
      cell[i][0] = update_rule(cell_old[i][0], sum_range(cell_old, i-1, i, 0, 1) - cell_old[i][0]);
      for ( j = 1; j < (n-1); j++)
      {
        cell[i][j] = update_rule(cell_old[i][j], sum_range(cell_old, i-1, i, j-1, j+1) - cell_old[i][j]);
      }
      cell[i][n-1] = update_rule(cell_old[i][n-1], sum_range(cell_old, i-1, i, n-2, n-1) - cell_old[i][n-1]);
    }
    else
    {
      for ( j = 0; j < n; j++)
      {
        cell[i][j] = update_rule(cell_old[i][j], sum_range(cell_old, i-1, i+1, j-1, j+1) - cell_old[i][j]);
      }
    }
  }
}

void print_cell(int ** cell, int m, int n)
{
  int i, j ;
  for ( i = 0; i < n; i++)
  {
    printf("--");
  }
  printf("--\n");
  for ( i = 0; i < m; i++)
  {
    printf("|");
    for ( j = 0; j < n; j++)
    {
      if(cell[i][j] == 1)
      {
        printf("o ");
      }
      else
      {
        printf("  ");
      }
    }
    printf("|\n");
  }
  for ( i = 0; i < n; i++)
  {
    printf("--");
  }
  printf("--\n");
}

void cell_automation(int m, int n, int gen_max, double * out_vec)
{
  int i, j, trial_count, gen_count;
  double t0, t1, t_end, t_avg;
  int ** cell = imatrix(0, m-1, 0, n-1);
  int ** cell_old = imatrix(0, m-1, 0, n-1);
  gsl_rng * T = gsl_rng_alloc(gsl_rng_taus);

  zeromint_init(cell, 0, m-1, 0, n-1);
  for ( i = ((m/2) - 6); i < ((m/2) + 6); i++)
  {
    for ( j = ((n/2) - 6); j < ((n/2) + 6); j++)
    {
      if (gsl_rng_uniform(T) < 0.75)
      {
        cell[i][j] = 1;
      }
    }
  }
  t_end = 0;
  for ( gen_count = 1; gen_count <= gen_max; gen_count++)
  {
    imatrix_cpy(cell, 0, m-1, 0, n-1, cell_old, 0, m-1, 0, n-1);
    t0 = omp_get_wtime();
    #pragma omp parallel for
        for(i=0; i<m; i++)
        {
          update_row(cell_old, i, m, n, cell);
        }
    t1 = omp_get_wtime();
    t_end += (t1-t0);
      // if (gen_count%25 == 0)
      // {
      //   printf("Gen %d cell: \n", gen_count);
      //   print_cell(cell, m, n);
      // }
  }
  // printf("Trial: %d final cell: \n", gen_count);
  // print_cell(cell, m, n);

  t_avg = (t_end)/( (double) gen_max );
  printf("num threads: %d, n: %d, mean gen calc time: %f\n", omp_get_max_threads(), n, t_avg);

  out_vec[0] = (double) omp_get_max_threads();
  out_vec[1] = (double) n;
  out_vec[2] = t_avg;

  free_imatrix(cell, 0, m-1, 0, n-1);
  free_imatrix(cell_old, 0, m-1, 0, n-1);
  gsl_rng_free(T);
}

uint64_t next_prime(int * status_vec, uint64_t last_prime)
{
  int stop = 0;
  uint64_t index = last_prime+2;
  uint64_t return_index;

  while(stop == 0)
  {
    if (status_vec[index] == 0)
    {
      return_index = index;
      stop = 1;
    }
    else
    {
      index = index + 2;
    }
  }
  return return_index;
}

uint64_t * find_primes(uint64_t max) // dummy simple Sieve of Eratosthenes since I'm not clever enough to implement the optimizations
{
  uint64_t i, j, prime_it;
  int stop = 0;
  // -1 denotes unsearched, 1 denotes prime, 0 denotes not prime
  int * status_vec = ivector(0, max);
  status_vec[0] = 0;
  status_vec[1] = 0;
  status_vec[2] = 1;
  status_vec[3] = 1;
  for ( i = 4; i <= max; i++) status_vec[i] = -1;

  prime_it = 3;
  while(stop == 0)
  {
    if (prime_it*prime_it > max)
    {
      stop = 1;
    }
    else
    {
      prime_it = next_prime(status_vec, prime_it);
      for ( i = 2*prime_it; i <= max; i = i + prime_it)
      {
        status_vec[i] = 0;
      }
    }
  }

  for ( i = 0; i < count; i++) ;



}
