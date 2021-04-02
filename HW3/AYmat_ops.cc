#include "AYmat_ops.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}
void AYmat_mul(AYmat * A1, AYmat * A2, AYmat * A3)
{
  if (A1->N != A2->M)
  {
    printf("Matrix mult failed: A1 N and A2 M not equal.\n");
  }
  else
  {
    if (A3->M != A1->M)
    {
      printf("Matrix mult failed: A1 M and A3 M not equal.\n");
    }
    else
    {
      if (A3->N != A2->N)
      {
        printf("Matrix mult failed: A2 N and A3 N not equal.\n");
      }
      else
      {
        dmatrix_mult(A1->A, 0, A1->M-1, 0, A1->N-1, A2->A, 0, A2->M-1, 0, A2->N-1, A3->A);
      }
    }
  }
}

void AYmat_mul_Strass(AYmat * A1, AYmat * A2, AYmat * A3)
{
  if (A1->N != A2->M)
  {
    printf("Matrix mult failed: A1 N and A2 M not equal.\n");
  }
  else
  {
    if (A3->M != A1->M)
    {
      printf("Matrix mult failed: A1 M and A3 M not equal.\n");
    }
    else
    {
      if (A3->N != A2->N)
      {
        printf("Matrix mult failed: A2 N and A3 N not equal.\n");
      }
      else
      {
        Strass_recurse(A1->A, A2->A, A3->A, A1->M);
      }
    }
  }
}

void Strass_recurse(double ** A, double ** B, double ** C, int N)
{
  if ((N > 2) && (N%2 == 0) )
  {
    int i, j;
    int n = N/2;
    int nn = n*n;

    double ** A00 = A;
    double ** A01 = (double**)malloc((n)*sizeof(double *));
    double ** A10 = A+N/2;
    double ** A11 = (double**)malloc((n)*sizeof(double *));

    double ** B00 = B;
    double ** B01 = (double**)malloc((n)*sizeof(double *));
    double ** B10 = B+N/2;
    double ** B11 = (double**)malloc((n)*sizeof(double *));

    double ** work1 = dmatrix(0, N/2-1, 0, N/2-1);

    for ( i = 0; i < n; i++)
    {
      A01[i] = A[i] + n;
      A11[i] = A10[i] + n;

      B01[i] = B[i] + n;
      B11[i] = B10[i] + n;
    }

    double ** Q1= dmatrix(0, n-1, 0, n-1); // Q1 gen
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work1[i][j] = A10[i][j] + A11[i][j];} }
    Strass_recurse(work1, B00, Q1, n);

    double ** Q2= dmatrix(0, n-1, 0, n-1);
    // Q2 gen
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work1[i][j] = B01[i][j] - B11[i][j];} }
    Strass_recurse(A00, work1, Q2, n);

    double ** Q3= dmatrix(0, n-1, 0, n-1);
    // Q3 gen
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work1[i][j] = B10[i][j] - B00[i][j];} }
    Strass_recurse(A11, work1, Q3, n);

    double ** Q4= dmatrix(0, n-1, 0, n-1);
    // Q4 gen
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work1[i][j] = A00[i][j] + A01[i][j];} }
    Strass_recurse(work1, B11, Q4, n);

    double ** work2 = dmatrix(0, n-1, 0, n-1);

    double ** Q0= dmatrix(0, n-1, 0, n-1); // Q0 gen
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work1[i][j] = A00[i][j] + A11[i][j];} }
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work2[i][j] = B00[i][j] + B11[i][j];} }
    Strass_recurse(work1, work2, Q0, n);

    double ** Q5= dmatrix(0, n-1, 0, n-1);
    // Q5 gen
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work1[i][j] = A10[i][j] - A00[i][j];} }
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work2[i][j] = B00[i][j] + B01[i][j];} }
    Strass_recurse(work1, work2, Q5, n);

    double ** Q6= dmatrix(0, n-1, 0, n-1);
    // Q6 gen
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work1[i][j] = A01[i][j] - A11[i][j];} }
    for ( i = 0; i < n; i++){ for ( j = 0; j < n; j++){work2[i][j] = B10[i][j] + B11[i][j];} }
    Strass_recurse(work1, work2, Q6, n);

    for ( i = 0; i < n; i++)
    {
      for ( j = 0; j < n; j++)
      {
        C[i][j] = Q0[i][j] + Q3[i][j] - Q4[i][j] + Q6[i][j];
        C[i+n][j] = Q1[i][j] + Q3[i][j];
        C[i][j+n] = Q2[i][j] + Q4[i][j];
        C[i+n][j+n] = Q0[i][j] + Q2[i][j] - Q1[i][j] + Q5[i][j];
      }
    }

    free_dmatrix(work1, 0, n-1, 0, n-1);
    free_dmatrix(work2, 0, n-1, 0, n-1);
    free_dmatrix(Q0, 0, n-1, 0, n-1);
    free_dmatrix(Q1, 0, n-1, 0, n-1);
    free_dmatrix(Q2, 0, n-1, 0, n-1);
    free_dmatrix(Q3, 0, n-1, 0, n-1);
    free_dmatrix(Q4, 0, n-1, 0, n-1);
    free_dmatrix(Q5, 0, n-1, 0, n-1);
    free_dmatrix(Q6, 0, n-1, 0, n-1);

    free(A01);
    free(A11);
    free(B01);
    free(B11);

  }
  else
  {
    dmatrix_mult(A, 0, N-1, 0, N-1, B, 0, N-1, 0, N-1, C);
  }
}

void AYmat_mul_BLAS(AYmat * A1, AYmat * A2, AYmat * A3)
{
  A1->gen_transpose();
  A2->gen_transpose();
}
