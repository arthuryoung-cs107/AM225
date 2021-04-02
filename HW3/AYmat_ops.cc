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
