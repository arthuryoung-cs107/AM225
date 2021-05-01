#include <cstdlib>

#include "cubic_1d_alt.hh"
#include "blas.h"

cubic_1d_alt::cubic_1d_alt(int n_) : conj_grad(n_+5), n(n_+5), n_full(n_+5), x_full(new double[n]),
h(1.0/( (double) n_ )), q(new quadrat(7)), node_pos(new double[n_full]),
omega(new double[2]), omega_1(new double[2]),
a_vals(dmatrix(0, n_full-1, -3, 3)), a_ind(imatrix(0, n_full-1, -3, 3)),
bounds(dmatrix(0, n_full-1, -3, 3)), a_short(dmatrix(0, n-1, -3, 3)),
b_full(new double[n_full])
{
  int i, j;

  omega[0] = 1.0;
  omega[1] = 2.0;
  omega_1[0] = omega[0] - 2.0*h;
  omega_1[1] = omega[1] + 2.0*h;

  zerom_init(a_vals, 0, n_full-1, -3, 3);
  zerom_init(a_short, 0, n-1, -3, 3);
  zeromint_init(a_ind, 0, n_full-1, -3, 3);
  for ( i = 0; i < n_full; i++)
  {
    node_pos[i] = omega_1[0] + h*( (double) i );
    bounds[i][0] = max(omega_1[0], node_pos[i] - 2.0*h);
    bounds[i][1] = min(omega_1[1], node_pos[i] + 2.0*h);
  }
  for ( i = 0; i < n_full; i++)
  {
    for ( j = -3; j <= 3; j++)
    {
      if (i + j > -1 && i + j < n_full)
      {
        a_ind[i][j] = 1;
      }
    }
  }
  assemble_a();
}

cubic_1d_alt::~cubic_1d_alt()
{
  delete q;
  delete [] node_pos;
  delete [] omega;
  delete [] omega_1;
  delete [] x_full;
  delete [] b_full;

  free_dmatrix(a_vals, 0, n_full-1, -3, 3);
  free_dmatrix(a_short, 0, n-1, -3, 3);
  free_imatrix(a_ind, 0, n_full-1, -3, 3);
  free_dmatrix(bounds, 0, n_full-1, -3, 3);
}

/** Performs multiplication on a vector by the stiffness matrix. */
void cubic_1d_alt::mul_A(double *in,double *out)
{
  int i, k;
  double acc;

  for ( i = 0; i < n; i++)
  {
    acc = 0.0;
    for ( k = -3; k <= 3; k++)
    {
      if (a_ind[i][k] == 1)
      {
        acc += a_vals[i][k]*in[i+k];
      }
    }
    out[i] = acc;
  }

}

void cubic_1d_alt::assemble_b()
{
  int i, j, k;
  double acc, val, C1, C2, C3, C4;
  for ( i = 0; i < n_full; i++)
  {
    acc = 0.0;
    C1 = bounds[i][0];
    C2 = bounds[i][1];
    C3 = (C2-C1)/2.0;
    C4 = (C2+C1)/2.0;
    for ( j = 0; j < q->n; j++)
    {
      acc += q->w[j]*(phi_C1(C3*(q->x[j]) + C4, i))*(f_source(C3*(q->x[j]) + C4));
    }
    b[i] = acc*C3;
    b[i] += 2.0*g*(phi_C1(omega[1], i));
  }
}

void cubic_1d_alt::assemble_a()
{
  int i, j, k;
  double acc, val, C1, C2, C3, C4, phi0, phi1, phi2;

  for ( i = 0; i < n_full; i++)
  {
    for ( k = -3; k <= 3; k++)
    {
      if (a_ind[i][k] == 1)
      {
        acc = 0.0;
        C1 = max(bounds[i][0], bounds[i+k][0]);
        C2 = min(bounds[i][1], bounds[i+k][1]);
        C3 = (C2-C1)/2.0;
        C4 = (C2+C1)/2.0;
        for ( j = 0; j < q->n; j++)
        {
          val = grad_phi_C1(C3*q->x[j] + C4, i)*grad_phi_C1(C3*q->x[j] + C4, i+k);
          acc += q->w[j]*(val)*(C3*q->x[j] + C4);
        }
        a_vals[i][k] = acc*C3;
      }
    }
  }

  phi0 = phi_C1(omega[0], 1);
  phi1 = phi_C1(omega[0], 2);
  phi2 = phi_C1(omega[0], 3);

  a_vals[1][0] -= (1.0/6.0)*phi0;
  a_vals[1][1] -= (1.0/6.0)*phi1;
  a_vals[1][2] -= (1.0/6.0)*phi2;

  a_vals[2][-1] -= (2.0/3.0)*phi0;
  a_vals[2][0]  -= (2.0/3.0)*phi1;
  a_vals[2][1]  -= (2.0/3.0)*phi2;

  a_vals[3][-2] -= (1.0/6.0)*phi0;
  a_vals[3][-1] -= (1.0/6.0)*phi1;
  a_vals[3][0]  -= (1.0/6.0)*phi2;

}

double cubic_1d_alt::phi_C1(double x_in, int i)
{
  double x_out=0.0;

  if (x_in < bounds[i][0] || x_in > bounds[i][1] )
  {
    return x_out;
  }
  else
  {
    double x_eval = (abs(x_in - node_pos[i]))/(h);
    if ( x_eval < 1.0 )
    {
      double del = 1.0 - x_eval;
      x_out = 0.25*(1.0 + 3.0*(del) + 3.0*(del*del) - 3.0*(del*del*del) );
    }
    else
    {
      double del = 2.0 - x_eval;
      x_out = 0.25*( del*del*del );
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
  if (x_in < bounds[i][0] || x_in > bounds[i][1] )
  {
    return x_out;
  }
  else
  {
    double x_eval = (abs(x_in - node_pos[i]))/(h);
    if ( x_eval < 1.0 )
    {
      double del = 1.0 - x_eval;
      x_out = -0.75 - (0.25)*6.0*del + (0.25)*9.0*del*del;
    }
    else
    {
      double del = 2.0 - x_eval;
      x_out = -0.75*( del*del );
    }
    if (x_in < node_pos[i])
    {
      x_out *= -1.0;
    }
    x_out*= (1.0/h);
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

  for ( i = 0; i < N_test; i++)
  {
    sol_out[i][0] = test_coords[i];
    acc = 0.0;
    for ( j = 0; j < n; j++)
    {
      if (test_coords[i] > bounds[j][0] && test_coords[i] < bounds[j][1])
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
