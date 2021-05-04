#include <cstdlib>

#include "cubic_1d_alt_C1.hh"
#include "blas.h"

cubic_1d_alt_C1::cubic_1d_alt_C1(int n_) : conj_grad(2*(n_)+1),
n_in(n_), n(2*(n_)+1), h(1.0/((double)n_)),
q(new quadrat(7)), q2(new quadrat(15)), node_pos(new double[n_+1]), omega(new double[2]), a_vals(dmatrix(0, n-1, -2, 2)),
a_ind(imatrix(0, n-1, -2, 2)), bounds(dmatrix(0, n-1, 0, 1))
{
  int i, j;

  omega[0] = 1.0;
  omega[1] = 2.0;

  zerom_init(a_vals, 0, n-1, -2, 2);
  zeromint_init(a_ind, 0, n-1, -2, 2);
  for ( i = 0; i <= n_in; i++)
  {
    node_pos[i] = omega[0] + h*( (double) i );
  }

  bounds[0][0] = omega[0]; bounds[0][1] = omega[0] + h;
  j = 1;
  for ( i = 1; i < n_in; i++)
  {
    bounds[j][0] = bounds[j+1][0] = node_pos[i]-h;
    bounds[j][1] = bounds[j+1][1] = node_pos[i]+h;
    j+=2;
  }
  bounds[n-2][0] = bounds[n-1][0] = omega[1]-h;
  bounds[n-2][1] = bounds[n-1][1] = omega[1];

  a_ind[0][0] = a_ind[0][1] = a_ind[0][2] = 1;
  a_ind[1][-1] = a_ind[1][0] = a_ind[1][1] = a_ind[1][2] = 1;
  for ( i = 2; i < n-2; i++)
  {
    for ( j = -2; j <= 2; j++)
    {
      a_ind[i][j] = 1;
    }
  }
  a_ind[n-2][1] = a_ind[n-2][0] = a_ind[n-2][-1] = a_ind[n-2][-2] = 1;
  a_ind[n-1][0] = a_ind[n-1][-1] = a_ind[n-1][-2] = 1;

  assemble_a();

}

cubic_1d_alt_C1::~cubic_1d_alt_C1()
{
  delete q;
  delete q2;
  delete [] node_pos;
  delete [] omega;

  free_dmatrix(a_vals, 0, n-1, -2, 2);
  free_imatrix(a_ind, 0, n-1, -2, 2);

}

/** Performs multiplication on a vector by the stiffness matrix. */
void cubic_1d_alt_C1::mul_A(double *in,double *out)
{
  int i, k;
  double acc;

  for ( i = 0; i < n; i++)
  {
    acc = 0.0;
    for ( k = -2; k <= 2; k++)
    {
      if (a_ind[i][k] == 1)
      {
        acc += a_vals[i][k]*in[i+k];
      }
    }
    out[i] = acc;
  }

}

void cubic_1d_alt_C1::assemble_b()
{
  int i, j, k;
  double acc, val, C1, C2, C3, C4;
  for ( i = 0; i < n; i++)
  {
    acc = 0.0;
    C1 = bounds[i][0];
    C2 = bounds[i][1];
    C3 = (C2-C1)/2.0;
    C4 = (C2+C1)/2.0;
    for ( j = 0; j < q2->n; j++)
    {
      acc += q2->w[j]*(phi_C1(C3*q2->x[j] + C4, i))*(f_source(C3*q2->x[j] + C4));
    }
    b[i] = acc*C3;
    b[i] += 2.0*g*(phi_C1(omega[1], i));
  }

}

void cubic_1d_alt_C1::assemble_a()
{
  int i, j, k;
  double acc, val, C1, C2, C3, C4, phi0, phi1, phi2;

  for ( i = 0; i < n; i++)
  {
    for ( k = -2; k <= 2; k++)
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

}

double cubic_1d_alt_C1::phi_C1(double x_in, int i)
{
  double x_out=0.0;

  if (x_in < bounds[i][0] || x_in > bounds[i][1] )
  {
    return x_out;
  }
  else
  {
    if (i > 0 && i < n-2 ) // is interior function?
    {
      if ( i%2 == 1) // is phi_odd?
      {
        int j = (i+1)/2;
        double x_eval = abs(x_in - node_pos[j])/h;
        x_out = (x_eval - 1.0)*(x_eval - 1.0)*x_eval; // phi_1, phi_3
        if (x_in < node_pos[j]) {x_out*=-1.0 ;}
      }
      else // is phi_even.
      {
        int j = (i+1)/2;
        double x_eval = abs(x_in - node_pos[j])/h;
        x_out = 2.0*(x_eval*x_eval*x_eval)-3.0*(x_eval*x_eval)+1.0; // phi_2, phi_0
      }
      return x_out;
    }
    else // is boundary function.
    {
      if (i == 0)
      {
        double x_eval = abs(x_in - node_pos[0])/h;
        x_out = (x_eval - 1.0)*(x_eval - 1.0)*x_eval; // phi_1
      }
      else
      {
        if (i == n-2)
        {
          double x_eval = abs(x_in - node_pos[n_in-1])/h;
          x_out = -(x_eval*x_eval)*(1.0 - x_eval); // phi_3
        }
        else
        {
          if (i == n-1)
          {
            double x_eval = abs(x_in - node_pos[n_in-1])/h;
            x_out = (x_eval*x_eval)*(3.0 - 2.0*x_eval); // phi_2
          }
          else
          {
            printf("we got a problem\n");
          }
        }
      }
      return x_out;
    }
  }
}

double cubic_1d_alt_C1::f_source(double xx)
{
  const double o=5.0*M_PI;
  return (-exp(1.0-xx)*(o*(1.0-2.0*xx)*cos(o*xx)+((1.0-o*o)*xx-1.0)*sin(o*xx)));
}

double cubic_1d_alt_C1::grad_phi_C1(double x_in, int i)
{
  double x_out=0.0;

  if (x_in < bounds[i][0] || x_in > bounds[i][1] )
  {
    return x_out;
  }
  else
  {
    if (i > 0 && i < n-2 ) // is interior function?
    {
      if ( i%2 == 1) // is phi_odd?
      {
        int j = (i+1)/2;
        double x_eval = abs(x_in - node_pos[j])/h;
        x_out = (x_eval - 1.0)*(x_eval - 1.0) + 2.0*(x_eval)*(x_eval - 1.0); // phi_1, phi_3
      }
      else // is phi_even.
      {
        int j = (i+1)/2;
        double x_eval = abs(x_in - node_pos[j])/h;
        x_out = 6.0*(x_eval*x_eval) - 6.0*x_eval; // phi_2, phi_0
        if (x_in < node_pos[j]) {x_out*=-1.0 ;}
      }
      x_out*= (1.0/h);
      return x_out;
    }
    else // is boundary function.
    {
      if (i == 0)
      {
        double x_eval = abs(x_in - node_pos[0])/h;
        x_out = (x_eval - 1.0)*(x_eval - 1.0) + 2.0*(x_eval)*(x_eval - 1.0); // phi_1
      }
      else
      {
        if (i == n-2)
        {
          double x_eval = abs(x_in - node_pos[n_in-1])/h;
          x_out = -2.0*x_eval*(1.0-x_eval) + (x_eval*x_eval);  // phi_3
        }
        else
        {
          if (i == n-1)
          {
            double x_eval = abs(x_in - node_pos[n_in-1])/h;
            x_out = 2.0*(x_eval)*(3.0 - 2.0*x_eval)-2.0*(x_eval*x_eval); //phi_2
          }
          else
          {
            printf("we got a problem\n");
          }
        }
      }
      x_out*= (1.0/h);
      return x_out;
    }
  }
}

void cubic_1d_alt_C1::write_out(char prefix[], int N_test)
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
      acc += x[j]*phi_C1(test_coords[i], j);
    }
    sol_out[i][1] = acc;
  }
  fprintf_matrix(sol_out, N_test, 2, prefix);

  delete [] test_coords;
  free_dmatrix(sol_out, 0, N_test-1, 0, 1);
}
