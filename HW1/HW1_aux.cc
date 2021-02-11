#include "HW1_aux.hh"
#include "omp.h"

extern "C" {
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

/** Custom random number generator based of "Ran" routine in Numerical Recipes
 * by Press et al. */
class custom_rng {
    public:
        unsigned long a,b,c;
        custom_rng(unsigned long seed) : b(4101842887655102017L), c(1) {
            if(sizeof(unsigned long)<8) {
                fputs("Error: 'unsigned long' type too short\n",stderr);
                exit(1);
            }
            a=seed^b;int64();
            b=a;int64();
            c=b;int64();
        }
        unsigned long int64() {
            a=a*2862933555777941757L+7046029254386353087L;
            b^=b>>17;b^=b<<31;b^=b>>8;
            c=4294957665U*(c&0xffffffff)+(c>>32);
            unsigned long d=a^(a<<21);
            d^=d>>35;d^=d<<4;
            return (d+b)^c;
        }
        inline double doub() {
            return 5.42101086242752217E-20*int64();
        }
};

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

void cell_automaton(int m, int n, int number_threads)
{
  int i, j, gen_count;
  int gen_max = 150;
  uint64_t start_seed = 123456789;
  uint64_t seed = lcg_fwd(start_seed, 100);

  int ** cell = imatrix(0, m-1, 0, n-1);
  int ** cell_old = imatrix(0, m-1, 0, n-1);

  zeromint_init(cell, 0, m-1, 0, n-1);

  for ( i = ((m/2) - 6); i < ((m/2) + 6); i++)
  {
    for ( j = ((n/2) - 6); j < ((n/2) + 6); j++)
    {
      if (random_uni(0, 1, &seed) < 0.75)
      {
        cell[i][j] = 1;
      }
    }
  }

  for ( gen_count = 0; gen_count < gen_max; gen_count++)
  {
    imatrix_cpy(cell, 0, m-1, 0, n-1, cell_old , 0, m-1, 0, n-1);
    #pragma omp parallel for num_threads(number_threads)
        for(i=0; i<m; i++)
        {
          update_row(cell_old, i, m, n, cell);
        }
      if (gen_count%25 == 0)
      {
        printf("Gen %d cell: \n", gen_count);
        print_cell(cell, m, n);
      }
  }

}
