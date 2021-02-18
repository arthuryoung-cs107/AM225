#include "HW1_aux.hh"
#include "omp.h"


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

long next_prime(short * status_vec, long last_prime)
{
  int stop = 0;
  long index = last_prime+2;
  long return_index;

  while(stop == 0)
  {
    if (status_vec[index] == -1)
    {
      status_vec[index] = 1; // is prime
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

Primes * find_primes(long check_prime)
{
  struct Prime_struct * return_primes = (Primes *)malloc(sizeof(Primes));
  long i, j, prime_it, num_primes;
  int stop = 0;
  // -1 denotes unsearched, 1 denotes prime, 0 denotes not prime
  short * status_vec = (short *)malloc((check_prime + 1)*sizeof(short));
  status_vec[0] = 0;
  status_vec[1] = 0;
  status_vec[2] = 1;
  for ( i = 3; i <= check_prime; i = i + 2) status_vec[i] = -1;
  for ( i = 4; i <= check_prime; i = i + 2) status_vec[i] = 0;

  prime_it = 1;
  num_primes = 1;
  while(stop == 0)
  {
    if (prime_it*prime_it > check_prime)
    {
      stop = 1;
    }
    else
    {
      num_primes++;
      prime_it = next_prime(status_vec, prime_it);
      for ( i = 2*prime_it; i <= check_prime; i = i + prime_it)
      {
        status_vec[i] = 0;
      }
    }
  }

  for ( i = prime_it+2; i <= check_prime; i = i + 2)
  {
    if (status_vec[i] == -1)
    {
      status_vec[i] = 1;
      num_primes++;
    }
  }

  long * prime_vec = (long *)malloc((num_primes)*sizeof(long));

  prime_vec[0] = 2;
  j = 1;
  for ( i = 3; i <= check_prime; i = i + 2)
  {
    if (status_vec[i] == 1)
    {
      prime_vec[j] = i;
      j++;
    }
  }

  return_primes->primes = prime_vec;
  return_primes->max = check_prime;
  return_primes->len = num_primes;

  free(status_vec);
  return return_primes;
}

Mersenne * Mersenne_expand(long n, long m_base)
{
  double pow_d;
  long i, local_pow, old_pow;
  long base = (long) pow((double) 2, (double) m_base);
  long base_1m = base - 1;
  long k = n/m_base; // integer division
  struct Mersenne_Struct * return_mersenne = (Mersenne *)malloc(sizeof(Mersenne));

  long * coeffs = (long * )malloc((k+1)*sizeof(long));
  return_mersenne->m = m_base;
  return_mersenne->k_max = k;
  return_mersenne->coeffs = coeffs;
  return_mersenne->base = base;


  old_pow = n;
  local_pow = old_pow - (m_base*k);
  pow_d = pow((double) 2, (double) local_pow);
  coeffs[k] = (long) pow_d - 1;

  for ( i = 0; i < k; i++) coeffs[i] = base_1m;

  return return_mersenne;
}

long general_divide(Mersenne * num_in, long divide)
{

  long i, whole, check1, check2;
  long acc = 0;
  long rem = 0;
  long base = (long) pow((double) 2, (double) num_in->m);
  for ( i = 0; i <= num_in->k_max; i++)
  {
    rem += ( (long) pow( (double) base, (double) i ) )*((num_in->coeffs[i])%divide);
    num_in->coeffs[i] = (num_in->coeffs[i])/divide;
    acc += rem/divide;
    rem = rem%divide;
  }

  num_in->coeffs[0] += acc;

  for ( i = 0; i <= num_in->k_max; i++)
  {
    if (num_in->coeffs[i] >= base)
    {
      num_in->coeffs[i+1]++;
      num_in->coeffs[i] -= base;
      i = i - 1;
    }
  }

  return rem;
}

void free_Primes(Primes * prime_structure)
{
  free(prime_structure->primes);
  free(prime_structure);
}

void free_Mersenne(Mersenne * mers_struct)
{
  free(mers_struct->coeffs);
  free(mers_struct);
}
