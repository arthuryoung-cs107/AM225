#include <cstdlib>

#include "cubic_1d_alt.hh"
#include "blas.h"

cubic_1d_alt::cubic_1d_alt(int n_) : conj_grad(n_), n(n_),
h(1.0/( (double) n )), q(new quadrat(7)), node_pos(new double[n]), omega(new double[2]),
a_vals(dmatrix(0, n-1, 0, 4)), a_count(new int[n])
{
  int i, j;

  omega[0] = 1.0;
  omega[1] = 2.0;
  zerom_init(a_vals, 0, n-1, 0, 2);

  for ( i = 0; i < n; i++)
  {
    node_pos[i] = omega[0] + h*( (double) i + 1 );
    a_count[i] = 0;
  }
  assemble_a();
}

cubic_1d_alt::~cubic_1d_alt()
{
  delete q;
  delete [] node_pos;
  delete [] omega;
  delete [] a_count;
  free_dmatrix(a_vals, 0, n-1, 0, 2);
}

/** Performs multiplication on a vector by the stiffness matrix. */
void cubic_1d_alt::mul_A(double *in,double *out)
{
  int i;
  out[0] = in[0]*(a_vals[0][0]) + in[1]*(a_vals[0][1]);
  for ( i = 1; i < n-1; i++)
  {
    out[i] = in[i-1]*(a_vals[i][0]) + in[i]*(a_vals[i][1]) + in[i+1]*(a_vals[i][2]);
  }
  out[n-1] = in[n-2]*(a_vals[n-1][0]) + in[n-1]*(a_vals[n-1][1]);
}

void cubic_1d_alt::assemble_b()
{
  int i, j;
  double acc, val, C1, C2, C3, C4;
  for ( i = 0; i < n-1; i++)
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

  acc = 0.0;
  C1 = node_pos[n-2];
  C2 = node_pos[n-1];
  C3 = (C2-C1)/2.0;
  C4 = (C2+C1)/2.0;
  for ( j = 0; j < q->n; j++)
  {
    acc += q->w[j]*(phi_C1(C3*q->x[j] + C4, n-1))*(f_source(C3*q->x[j] + C4));
  }
  b[n-1] = acc*C3;

  b[n-1]+=2.0*g;
}

void cubic_1d_alt::assemble_a()
{
  int i, j;
  double acc, val, C1, C2, C3, C4;

  acc = 0.0;
  C1 = node_pos[0]-h;
  C2 = node_pos[0]+h;
  C3 = (C2-C1)/2.0;
  C4 = (C2+C1)/2.0;
  for ( j = 0; j < q->n; j++)
  {
    val = grad_phi_C1(C3*q->x[j] + C4, 0);
    acc += q->w[j]*(val*val)*(C3*q->x[j] + C4);
    // acc += q->w[j]*(val*val);
  }
  a_vals[0][0] = acc*C3;

  acc = 0.0;
  C1 = node_pos[0];
  C2 = node_pos[0]+h;
  C3 = (C2-C1)/2.0;
  C4 = (C2+C1)/2.0;
  for ( j = 0; j < q->n; j++)
  {
    val = grad_phi_C1(C3*q->x[j] + C4, 0)*grad_phi_C1(C3*q->x[j] + C4, 1);
    acc += q->w[j]*(val)*(C3*q->x[j] + C4);
    // acc += q->w[j]*(val);
  }
  a_vals[0][1] = acc*C3;

  for ( i = 1; i < n-1; i++)
  {
    acc = 0.0;
    C1 = node_pos[i-1];
    C2 = node_pos[i];
    C3 = (C2-C1)/2.0;
    C4 = (C2+C1)/2.0;
    for ( j = 0; j < q->n; j++)
    {
      val = grad_phi_C1(C3*q->x[j] + C4, i)*grad_phi_C1(C3*q->x[j] + C4, i-1);
      acc += q->w[j]*(val)*(C3*q->x[j] + C4);
      // acc += q->w[j]*(val);
    }
    a_vals[i][0] = acc*C3;

    acc = 0.0;
    C1 = node_pos[i-1];
    C2 = node_pos[i+1];
    C3 = (C2-C1)/2.0;
    C4 = (C2+C1)/2.0;
    for ( j = 0; j < q->n; j++)
    {
      val = grad_phi_C1(C3*q->x[j] + C4, i);
      acc += q->w[j]*(val*val)*(C3*q->x[j] + C4);
      // acc += q->w[j]*(val*val);
    }
    a_vals[i][1] = acc*C3;

    acc = 0.0;
    C1 = node_pos[i];
    C2 = node_pos[i+1];
    C3 = (C2-C1)/2.0;
    C4 = (C2+C1)/2.0;
    for ( j = 0; j < q->n; j++)
    {
      val = grad_phi_C1(C3*q->x[j] + C4, i)*grad_phi_C1(C3*q->x[j] + C4, i+1);
      acc += q->w[j]*(val)*(C3*q->x[j] + C4);
      // acc += q->w[j]*(val);
    }
    a_vals[i][2] = acc*C3;
  }

  acc = 0.0;
  C1 = node_pos[n-2];
  C2 = node_pos[n-1];
  C3 = (C2-C1)/2.0;
  C4 = (C2+C1)/2.0;
  for ( j = 0; j < q->n; j++)
  {
    val = grad_phi_C1(C3*q->x[j] + C4, n-2)*grad_phi_C1(C3*q->x[j] + C4, n-1);
    acc += q->w[j]*(val)*(C3*q->x[j] + C4);
    // acc += q->w[j]*(val);
  }
  a_vals[n-1][0] = acc*C3;

  acc = 0.0;
  C1 = node_pos[n-2];
  C2 = node_pos[n-1];
  C3 = (C2-C1)/2.0;
  C4 = (C2+C1)/2.0;
  for ( j = 0; j < q->n; j++)
  {
    val = grad_phi_C1(C3*q->x[j] + C4, n-1);
    acc += q->w[j]*(val*val)*(C3*q->x[j] + C4);
    // acc += q->w[j]*(val*val);
  }
  a_vals[n-1][1] = acc*C3;

}

double cubic_1d_alt::phi_C1(double x_in, int i)
{
  double x_out=0.0;
  double x_eval = (abs(x_in - node_pos[i]))/h;

  if (i == n-1)
  {
    if (x_in > omega[1] || x_eval > 1.0)
    {
      return x_out;
    }
    else
    {
      x_out = 1.0  - x_eval - (x_eval*x_eval) + (x_eval*x_eval*x_eval);
      // x_out = 1.0 - 3.0*(x_eval*x_eval) + 2.0*(x_eval*x_eval*x_eval);
    }
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
  const double o=5.0*M_PI;
  return (-exp(1.0-xx)*(o*(1.0-2.0*xx)*cos(o*xx)+((1.0-o*o)*xx-1.0)*sin(o*xx)));
}

double cubic_1d_alt::grad_phi_C1(double x_in, int i)
{
  double x_out=0.0;
  double x_eval = (abs(x_in - node_pos[i]))/h;

  if (i == n-1)
  {
    if (x_in > omega[1] || x_eval > 1.0)
    {
      return x_out;
    }
    else
    {
      x_out = -(-1.0 - 2.0*(x_eval) + 3.0*(x_eval*x_eval))*(1.0/h);
      // x_out = -6.0*(-(x_eval) + (x_eval*x_eval) )*(1.0/h) ;
    }
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

void cubic_1d_alt::write_out(char prefix[], int N_test)
{
  int i, j;
  double acc;
  double del = (omega[1] - omega[0])/( (double) N_test - 1 );
  double ** sol_out = dmatrix(0, N_test-1, 0, 1);
  double * test_coords = new double[N_test];

  for ( i = 0; i < N_test; i++) test_coords[i] = omega[0] + ((double) i)*del;

  sol_out[0][0] = test_coords[0];
  sol_out[0][1] = 0.0;

  for ( i = 1; i < N_test; i++)
  {
    sol_out[i][0] = test_coords[i];
    acc = 0;
    for ( j = 0; j < n; j++)
    {
      if ( abs(test_coords[i] - node_pos[j]) < h )
      {
        acc += x[j]*phi_C1(test_coords[i], j);
      }
    }
    sol_out[i][1] = acc;
  }
  fprintf_matrix(sol_out, N_test, 2, prefix);

  delete [] test_coords;
  free_dmatrix(sol_out, 0, N_test-1, 0, 1);

}
