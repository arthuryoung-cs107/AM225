#include <cstdlib>

#include "cubic_1d_alt.hh"
#include "blas.h"

cubic_1d_alt::cubic_1d_alt(int n_) : conj_grad(n_), n(n_),
h(1.0/( (double) n )), q(new quadrat(10)), node_pos(new double[n]), omega(new double[2])
{
  double acc, val, val2, C1, C2, C3, C4;

  omega[0] = 1.0;
  omega[1] = 2.0;
  for (int i = 0; i < n; i++) node_pos[i] = omega[0] + h*( (double) i + 1 );

  acc = 0.0;
  C1 = node_pos[0]-h;
  C2 = node_pos[0]+h;
  C3 = (C2-C1)/2.0;
  C4 = (C2+C1)/2.0;
  for ( int j = 0; j < q->n; j++)
  {
    val = grad_phi_C1(C3*q->x[j] + C4, 0);
    acc += q->w[j]*(val*val);
  }
  a_C = acc*C3;

  acc = 0.0;
  C1 = node_pos[0]-h;
  C2 = node_pos[0]+h;
  C3 = (C2-C1)/2.0;
  C4 = (C2+C1)/2.0;
  for ( int j = 0; j < q->n; j++)
  {
    val = grad_phi_C1(C3*q->x[j] + C4, 0);
    val2 = grad_phi_C1(C3*q->x[j] + C4, 1);
    acc += q->w[j]*(val*val2);
  }
  a_C = acc*C3;
}

cubic_1d_alt::~cubic_1d_alt()
{
  delete q;
}

/** Performs multiplication on a vector by the stiffness matrix. */
void cubic_1d_alt::mul_A(double *in,double *out)
{

}

/** Computes the source vector in the linear system, which is based
 * on multiplying the source function by the mass matrix. */
void cubic_1d_alt::assemble_b()
{
  int i, j;
  double acc, x_loc, C1, C2, C3, C4;
  for ( i = 0; i < n; i++)
  {
    acc = 0.0;
    C1 = node_pos[i]-h;
    C2 = node_pos[i]+h;
    C3 = (C2-C1)/2.0;
    C4 = (C2+C1)/2.0;
    for ( j = 0; j < q->n; j++)
    {
      acc += q->w[j]*(phi_C1(C3*q->x[j] + C4, i))*(f_source(C3*q->x[j] + C4));
    }
    b[i] = acc*C3;
  }
}

double cubic_1d_alt::phi_C1(double x_in, int i)
{
  double x_out=0.0;
  double x_eval = (abs(x_in - node_pos[i]))/h;

  if (i == n-1 && x_in>2.0)
  {
    return x_out;
  }
  else
  {
    if (x_eval < 1.0)
    {
      x_out = 1.0 - 3.0*(x_eval*x_eval) + 2.0*(x_eval*x_eval*x_eval);
    }
    return x_out;
  }
}

double cubic_1d_alt::f_source(double xx)
{
  const double o=5*M_PI;
  return (-exp(1-xx)*(o*(1-2*xx)*cos(o*xx)+((1-o*o)*xx-1)*sin(o*xx)));
}

double cubic_1d_alt::grad_phi_C1(double x_in, int i)
{
  double x_out=0.0;
  double x_eval = (abs(x_in - node_pos[i]))/h;

  if (i == n-1 && x_in>2.0)
  {
    return x_out;
  }
  else
  {
    if (x_eval < 1.0)
    {
      x_out = 6.0*(-(x_eval) + (x_eval*x_eval) )*(1.0/h) ;
      if (x_in < node_pos[i])
      {
        x_out *= -1.0;
      }
    }
    return x_out;
  }
}
