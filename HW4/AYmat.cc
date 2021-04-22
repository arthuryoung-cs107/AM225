#include <cstdio>

#include "AYmat.hh"
#include "blas.h"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

// note: this is designed to play nicely with BLAS, hence we store column major
AYmat::AYmat(int M_, int N_): M(M_), N(N_), AT(dmatrix(0, N-1, 0, M-1))
{
  A_ptr = AT[0];
}

AYmat::~AYmat()
{
  free_dmatrix(AT, 0, N-1, 0, M-1);
}

void AYmat::print_mat()
{
  int i, j;
  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      printf("%f ", AT[j][i]);
    }
    printf("\n");
  }
}

void AYmat::init_123()
{
  int i, j, count;
  count = 1;
  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      AT[j][i] = (double) count;
      count++;
    }
  }
}

void AYmat::init_0()
{
  int i;
  for ( i = 0; i < N*M; i++) AT[0][i] = 0.0;
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
    AT[0][i] = check;
  }
}

void AYmat::set(int i, int j, double val) {AT[j][i] = val;}

double AYmat::get(int i, int j) {return AT[j][i];}

// take in two matrices, multiply A and B, set contents of this one, C
void AYmat::mult_set(AYmat * A_in, AYmat * B_in, double alpha, double beta )
{
  char trans = 'n';
  int inc = 1;
  if (B_in->N == 1) // matrix vector multiplication
  {
    dgemv_(&trans, &(A_in->M), &(A_in->N), &alpha, (A_in->A_ptr), &(A_in->M), (B_in->A_ptr), &inc, &beta, A_ptr, &inc);
  }
  else // matrix-matrix multiplication
  {
    dgemm_(&trans,&trans,&(A_in->M),&(N),&(B_in->M),&alpha,(A_in->A_ptr),&(A_in->M), (B_in->A_ptr),&(B_in->M),&beta,A_ptr,&(M));
  }
}

// take in two matrices, multiply this one and the A, set C with outcome
void AYmat::mult_put(AYmat * B_in, AYmat * C_in, double alpha, double beta )
{
  char trans = 'n';
  int inc = 1;
  if (B_in->N == 1) // matrix vector multiplication
  {
    dgemv_(&trans, &(M), &(N), &alpha, (A_ptr), &(M), (B_in->A_ptr), &inc, &beta, (C_in->A_ptr), &inc);
  }
  else // matrix-matrix multiplication
  {
    dgemm_(&trans,&trans,&(M),&(C_in->N),&(B_in->M),&alpha,(A_ptr),&(M), (B_in->A_ptr),&(B_in->M),&beta,(C_in->A_ptr),&(C_in->M));
  }
}
