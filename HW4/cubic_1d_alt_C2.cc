#include <cstdlib>

#include "cubic_1d_alt_C2.hh"
#include "blas.h"

cubic_1d_alt_C2::cubic_1d_alt_C2(int n_) : conj_grad(n_+3), n_in(n_), n(n_+3), n_full(n_+5), n_sol(n_+4),
x_sol(new double[n_sol]),
h(1.0/( (double) n_ )), q(new quadrat(7)), node_pos(new double[n_full]),
omega(new double[2]), omega_1(new double[2]), omega_2(new double[2]),
a_vals(dmatrix(0, n_full-1, -3, 3)), a_ind(imatrix(0, n_full-1, -3, 3)),
bounds(dmatrix(0, n_full-1, -3, 3)), a_short(dmatrix(0, n-1, -3, 3)),
b_full(new double[n_full]), ind(ivector(-1, n_full-1))
{
  int i, j;

  omega[0] = 1.0;
  omega[1] = 2.0;
  omega_1[0] = omega[0] - 1.0*h;
  omega_1[1] = omega[1] + 1.0*h;
  omega_2[0] = omega[0] - 2.0*h;
  omega_2[1] = omega[1] + 2.0*h;

  zerom_init(a_vals, 0, n_full-1, -3, 3);
  zerom_init(a_short, 0, n-1, -3, 3);
  zeromint_init(a_ind, 0, n_full-1, -3, 3);
  for ( i = 0; i < n_full; i++)
  {
    ind[i-1] = i;
    node_pos[i] = omega_2[0] + h*( (double) i );
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

cubic_1d_alt_C2::~cubic_1d_alt_C2()
{
  delete q;
  delete [] node_pos;
  delete [] omega;
  delete [] omega_1;
  delete [] omega_2;
  delete [] x_sol;
  delete [] b_full;

  free_dmatrix(a_vals, 0, n_full-1, -3, 3);
  free_dmatrix(a_short, 0, n-1, -3, 3);
  free_imatrix(a_ind, 0, n_full-1, -3, 3);
  free_ivector(ind, -1, n_in+3);
  free_dmatrix(bounds, 0, n_full-1, -3, 3);
}

/** Performs multiplication on a vector by the stiffness matrix. */
void cubic_1d_alt_C2::mul_A(double *in,double *out)
{
  int i, k;
  double acc;

  for ( i = 0; i < n; i++)
  {
    acc = 0.0;
    for ( k = -3; k <= 3; k++)
    {
      if (a_ind[i+2][k] == 1)
      {
        acc += a_short[i][k]*in[i+k];
      }
    }
    out[i] = acc;
  }

}

void cubic_1d_alt_C2::assemble_b()
{
  int i, j, k;
  double acc, val, C1, C2, C3, C4,
  a00, a01, a02, a03, am1_0, am1_1, am1_2, phi0, phi1, phi2,
  zeta01, zeta12, zeta02, zeta21, zeta0, zeta3;
  for ( i = 0; i < n_full; i++)
  {
    acc = 0.0;
    C1 = bounds[i][0];
    C2 = bounds[i][1];
    C3 = (C2-C1)/2.0;
    C4 = (C2+C1)/2.0;
    for ( j = 0; j < q->n; j++)
    {
      acc += q->w[j]*(phi_C2(C3*q->x[j] + C4, i))*(f_source(C3*q->x[j] + C4));
    }
    b_full[i] = acc*C3;
    b_full[i] += 2.0*g*(phi_C2(omega[1], i));
  }

  for ( i = 0; i < n; i++)
  {
    b[i] = b_full[i+2];
    for ( j = -3; j <= 3; j++)
    {
      a_short[i][j] = a_vals[i+2][j];
    }
  }


  // assert Dirichlet, condense system into n-2 dof

  a00 = a_vals[ind[0]][0];
  a01 = a_vals[ind[0]][1];
  a02 = a_vals[ind[0]][2];
  a03 = a_vals[ind[0]][3];

  am1_0 = a_vals[ind[-1]][1];
  am1_1 = a_vals[ind[-1]][2];
  am1_2 = a_vals[ind[-1]][3];

  phi0 = phi_C2(omega[0], ind[0]);
  phi1 = phi_C2(omega[0], ind[1]);
  phi2 = phi_C2(omega[0], ind[2]);

  zeta01 = am1_0 - (phi0*(am1_2/phi2));
  zeta12 = am1_1 - (phi1*(am1_2/phi2));

  zeta02 = am1_0 - (phi0*(am1_1/phi1));
  zeta21 = am1_2 - (phi2*(am1_1/phi1));

  zeta0 = a00 - (a01*zeta01/zeta12) - (a02*zeta02/zeta21);
  zeta3 = b_full[ind[0]] - (a01*b_full[ind[-1]]/zeta12) - (a02*b_full[ind[-1]]/zeta21);

  b[0] -= a01*b_full[ind[-1]]/zeta01;
  a_short[0][0] -= a01*zeta12/zeta01;

  b[1] -= a02*b_full[0]/zeta02;
  a_short[1][0] -= a02*zeta21/zeta02;

  b[2] -= a03*(zeta3)/zeta0;
  a_short[2][0] -= a03*a03/zeta0;

  a_ind[ind[1]][-3] = a_ind[ind[1]][-2] = a_ind[ind[1]][-1] = 0;
  a_ind[ind[2]][-3] = a_ind[ind[2]][-2] = 0;
  a_ind[ind[3]][-3] = 0;

  a_short[0][-3] = a_short[0][-2] = a_short[0][-1] = 0.0;
  a_short[1][-3] = a_short[1][-2] = 0.0;
  a_short[2][-3] = 0.0;

  // for ( i = 0; i < n; i++)
  // {
  //   printf("%d (%f) : ", i+1, node_pos[ind[i + 1]]);
  //   for ( k = 0; k < n; k++)
  //   {
  //     if (abs(i-k) < 4 )
  //     {
  //       printf("%f ", a_short[i][k-i]);
  //     }
  //     else
  //     {
  //       printf("%f ", 0.0);
  //     }
  //   }
  //   printf("|| %f \n", b[i]);
  // }
  //
  // getchar();
}

void cubic_1d_alt_C2::assemble_a()
{
  int i, j, k;
  double acc, val, C1, C2, C3, C4;

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
          val = grad_phi_C2(C3*q->x[j] + C4, i)*grad_phi_C2(C3*q->x[j] + C4, i+k);
          acc += q->w[j]*(val)*(C3*q->x[j] + C4);
        }
        a_vals[i][k] = acc*C3;
      }
    }
  }
}

double cubic_1d_alt_C2::phi_C2(double x_in, int i)
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

double cubic_1d_alt_C2::f_source(double xx)
{
  const double o=5.0*M_PI;
  return (-exp(1.0-xx)*(o*(1.0-2.0*xx)*cos(o*xx)+((1.0-o*o)*xx-1.0)*sin(o*xx)));
}

double cubic_1d_alt_C2::grad_phi_C2(double x_in, int i)
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

void cubic_1d_alt_C2::write_out(char prefix[], int N_test)
{
  int i, j;
  double acc;
  double del = (omega[1] - omega[0])/( (double) N_test - 1 );
  double ** sol_out = dmatrix(0, N_test-1, 0, 1);
  double * test_coords = new double[N_test];

  for ( i = 0; i < N_test; i++) test_coords[i] = omega[0] + ((double) i)*del;

  for ( i = 0; i < n; i++) x_sol[i+1] = x[i];
  x_sol[0] = -(phi_C2(1.0, ind[1])*x[0] +  phi_C2(1.0, ind[2])*x[1])/(phi_C2(1.0, ind[0])); // still valid

  for ( i = 0; i < N_test; i++)
  {
    sol_out[i][0] = test_coords[i];
    acc = 0.0;
    for ( j = 0; j < n_sol; j++)
    {
      acc += x_sol[j]*phi_C2(test_coords[i], ind[j]);
    }
    sol_out[i][1] = acc;
  }
  fprintf_matrix(sol_out, N_test, 2, prefix);

  delete [] test_coords;
  free_dmatrix(sol_out, 0, N_test-1, 0, 1);

}
