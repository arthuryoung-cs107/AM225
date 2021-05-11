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
AYmat::AYmat(int M_, int N_): M(M_), N(N_), AT(dmatrix(0, N-1, 0, M-1)), GSL_flag(0)
{
  A_ptr = AT[0];
}

AYmat::~AYmat()
{
  free_dmatrix(AT, 0, N-1, 0, M-1);
  if (GSL_flag == 1)
  {
    gsl_matrix_free(A_gsl);
  }
}

void AYmat::GSL_init()
{
  int i, j;
  if (GSL_flag == 0)
  {
    A_gsl = gsl_matrix_alloc(M, N);
    GSL_flag = 1;
  }

  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      gsl_matrix_set(A_gsl, i, j, AT[j][i]);
    }
  }
}

void AYmat::GSL_send()
{
  int i, j;
  if (GSL_flag == 0)
  {
    printf("GSL matrix not allocated.");
  }
  else
  {
    for ( i = 0; i < M; i++)
    {
      for ( j = 0; j < N; j++)
      {
        AT[j][i] = gsl_matrix_get(A_gsl, i, j);
      }
    }
  }
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

AYmat * AYmat::copy_gen()
{
  AYmat * X_out = new AYmat(M, N);
  memcpy(X_out->A_ptr, A_ptr, M*N*sizeof(double));
  return X_out;
}

void AYmat::copy_set(AYmat * X_in) { memcpy(X_in->A_ptr, A_ptr, M*N*sizeof(double)); }

AYmat * AYmat::transpose_gen()
{
  int i, j;
  AYmat * X_out = new AYmat(N, M);
  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      X_out->A_ptr[i*N + j] = AT[j][i];
    }
  }
  return X_out;
}

void AYmat::transpose_set(AYmat * X_in)
{
  int i, j;
  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      X_in->A_ptr[i*N + j] = AT[j][i];
    }
  }
}

void AYmat::add(AYmat * B_in, AYmat * C_in, double alpha, double beta )
{
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      C_in->AT[j][i] = AT[j][i] + alpha*B_in->AT[j][i] + beta*C_in->AT[j][i];
    }
  }
}

// take in one matrix, post multiply this matrix with input, return pointer to new matrix
AYmat * AYmat::mult_gen(AYmat * B_in, double alpha)
{
  double beta = 0.0;
  char trans = 'n';
  int inc = 1;

  AYmat * C_out = new AYmat(M, B_in->N);

  if (B_in->N == 1) // matrix vector multiplication
  {
    dgemv_(&trans, &(M), &(N), &alpha, (A_ptr), &(M), (B_in->A_ptr), &inc, &beta, (C_out->A_ptr), &inc);
  }
  else // matrix-matrix multiplication
  {
    dgemm_(&trans,&trans,&(M),&(C_out->N),&(B_in->M),&alpha,(A_ptr),&(M), (B_in->A_ptr),&(B_in->M),&beta,(C_out->A_ptr),&(C_out->M));
  }
  return C_out;
}

// take in two matrices, post multiply this one with B, set contents of C_in
void AYmat::mult_set(AYmat * B_in, AYmat * C_in, double alpha, double beta )
{
  char trans = 'n';
  int inc = 1;
  if (B_in->N == 1) // matrix vector multiplication
  {
    dgemv_(&trans, &(M), &(N), &alpha, (A_ptr), &(M), (B_in->A_ptr), &inc, &beta, C_in->A_ptr, &inc);
  }
  else // matrix-matrix multiplication
  {
    dgemm_(&trans,&trans,&(M),&(C_in->N),&(B_in->M),&alpha,(A_ptr),&(M), (B_in->A_ptr),&(B_in->M),&beta,(C_in->A_ptr),&(C_in->M));
  }
}

double AYmat::inner(AYmat * B_in)
{
  int i, j;
  double sum = 0;
  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      sum += (AT[j][i] * B_in->AT[j][i]);
    }
  }
  return sum;
}

double AYmat::norm_frob()
{
  int i,j;
  double out = 0.0;

  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      out += (AT[j][i])*(AT[j][i]);
    }
  }
  return sqrt(out);
}

double AYmat::norm_1()
{
  int i, j;
  double out = 0.0;

  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      out += abs(AT[j][i]);
    }
  }

  return out;
}

void AYmat::svd(gsl_vector * S_in, gsl_matrix * V_in, gsl_vector * work) {GSL_init(); gsl_linalg_SV_decomp(A_gsl, V_in, S_in, work);}

AYmat * aysml_read(char name[])
{
  int i, j, M, N;
  double data;
  size_t success;
  char extract[1000];
  char aysml_specfile[200];
  char aydat_specfile[200];
  memset(aysml_specfile, 0, 199);
  memset(aydat_specfile, 0, 199);
  snprintf(aysml_specfile, 200, "%s.aysml", name);
  snprintf(aydat_specfile, 200, "%s.aydat", name);

  FILE * aysml_stream = fopen(aysml_specfile, "r");

  char * check = fgets(extract, sizeof(extract), aysml_stream);
  M = atoi(extract);
  int next = (int) log10((double) M) + 1;
  N = atoi(extract+next);

  AYmat * out = new AYmat(M, N);
  fclose(aysml_stream);

  FILE * aydat_stream = fopen(aydat_specfile, "r");

  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      success = fread(&data, sizeof(double), 1, aydat_stream );
      out->set(i, j, data);
    }
  }
  fclose(aydat_stream);

  return out;
}
