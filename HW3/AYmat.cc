#include <cstdio>

#include "AYmat.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

AYmat::AYmat(int M_, int N_): M(M_), N(N_), A(dmatrix(0, M-1, 0, N-1))
{}

AYmat::~AYmat()
{
  free_dmatrix(A, 0, M-1, 0, N-1);
}

void AYmat::print_mat()
{
  int i, j;
  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      printf("%f ", A[i][j]);
    }
    printf("\n");
  }
}

void AYmat::init_123()
{
  int i;
  for ( i = 0; i < N*M; i++) A[0][i] = (double) (i + 1);
}

void AYmat::init_randuni()
{
  int i;
  uint64_t jump = 1000;
  uint64_t seed  = (uint64_t) rand();
  uint64_t carry = lcg_fwd(seed, jump);
  double low = -1.0;
  double high = 1.0;

  for ( i = 0; i < N*M; i++)
  {
    double check = knuth_random_uni(low, high, &carry);
    A[0][i] = check;
  }
}

void AYmat::gen_transpose()
{
  int i, j;
  AT = dmatrix(0, N-1, 0, M-1);

  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      AT[j][i] = 
    }
  }

}
