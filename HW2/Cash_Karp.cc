#include <cstdio>
#include <cmath>

#include "Cash_Karp.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

Cash_Karp::Cash_Karp(int dof_in): dof(dof_in), evals(0), h_init(0.01), fac(0.9), facmax(3.0), facmin(1.0/3.0),
k1(new double[dof]), k2(new double[dof]), k3(new double[dof]), k4(new double[dof]), k5(new double[dof]), k6(new double[dof]),
k_old(new double[dof]), k_new(new double[dof]),
w_it(new double[dof]), y_init(new double[dof]), y0(new double[dof]), y1(new double[dof]), y2(new double[dof]), yw(new double[dof]), yhat1(new double[dof]), yhat2(new double[dof]),
Theta(gsl_matrix_alloc(4, 4)), a_full(gsl_matrix_alloc(dof, 6)), b_vec(gsl_vector_alloc(4)), a_vec(gsl_vector_alloc(4)), perm_it(gsl_permutation_alloc(4)) {}

Cash_Karp::~Cash_Karp()
{
  delete [] k1;
  delete [] k2;
  delete [] k3;
  delete [] k4;
  delete [] k5;
  delete [] k6;
  delete [] k_old;
  delete [] k_new;
  delete [] w_it;
  delete [] y_init;
  delete [] y0;
  delete [] y1;
  delete [] y2;
  delete [] yw;
  delete [] yhat1;
  delete [] yhat2;

  gsl_matrix_free(Theta);
  gsl_matrix_free(a_full);
  gsl_vector_free(b_vec);
  gsl_vector_free(a_vec);
  gsl_permutation_free(perm_it);
}

void Cash_Karp::step(double * y_in, double * y_k, double t)
{
  int i;
  const double c2 = 0.2,
               c3 = 0.3,
               c4 = 0.6,
               c5 = 1.0,
               c6 = 7.0/8.0,

               b1 = 37.0/378.0,
               b2 = 0.0,
               b3 = 250.0/621.0,
               b4 = 125.0/594.0,
               b5 = 0.0,
               b6 = 512.0/1771.0,

               a21 = 0.2,
               a31 = 0.075,
               a32 = 0.225,
               a41 = 0.3,
               a42 = -0.9,
               a43 = 1.2,
               a51 = -11.0/54.0,
               a52 = 2.5,
               a53 = -70.0/27.0,
               a54 = 35.0/27.0,
               a61 = 1631.0/55296.0,
               a62 = 175.0/512.0,
               a63 = 575.0/13824.0,
               a64 = 44275.0/110592.0,
               a65 = 253.0/4096.0;

  for (i = 0; i < dof; i++) w_it[i] = y_in[i] + h*a21*k1[i];
  eval(t + c2*h, w_it, k2);

  for (i = 0; i < dof; i++) w_it[i] = y_in[i] + h*(a31*k1[i] + a32*k2[i]);
  eval(t + c3*h, w_it, k3);

  for (i = 0; i < dof; i++) w_it[i] = y_in[i] + h*(a41*k1[i] + a42*k2[i] + a43*k3[i]);
  eval(t + c4*h, w_it, k4);

  for (i = 0; i < dof; i++) w_it[i] = y_in[i] + h*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
  eval(t + c5*h, w_it, k5);

  for (i = 0; i < dof; i++) w_it[i] = y_in[i] + h*(a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
  eval(t + c6*h, w_it, k6);

  // add em up
  for (i = 0; i < dof; i++) y_k[i] = y_in[i] + h*( b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i]);
}

void Cash_Karp::extrap(double * yw_in , double * y_in1, double * y_in2, double * y_out1, double * y_out2)
{
  for (int i = 0; i < dof; i++)
  {
    y_out1[i] = y_in1[i] + ( (y_in2[i] - yw_in[i])/(62.0) ) ;
    y_out2[i] = y_in2[i] + ( (y_in2[i] - yw_in[i])/(31.0) ) ;
  }
}

double Cash_Karp::h_select(double * y_in, double h_in)
{
  int i;
  double err, h_new, h_it;
  h_it = h_in;
  for ( i = 0; i < dof; i++) k_old[i] = k_new[i];
  do
  {
    for ( i = 0; i < dof; i++) k1[i] = k_old[i];
    h = h_it;

    step(y_in, yw, t_it);

    h = h_it/(2.0);
    step(y_in, y1, t_it);

    eval(t_it+h, y1, k1);
    step(y1, y2, t_it+h);

    extrap(yw, y1, y2, yhat1, yhat2);

    err = err_est(y_in, y2, yhat2);
    h_new = h_it*min_d(facmax, max_d( facmin, fac*pow(err, (-1.0)/(6.0))));

    if (std::isnan(err))
    {
      err = 2.0;
      h_it = h_it*facmin;
    }
    else
    {
      if (err > 1.0)
      {
        h_it = h_new ;
      }
    }
  } while(err > 1.0);
  t_it += h_it;
  return h_new;
}

double Cash_Karp::err_est(double * y0_in, double * y2_in, double * y2hat_in)
{
  double sc, diff;
  double err_it = 0;
  for (int i = 0; i < dof; i++)
  {
    sc = a_tol + r_tol*max_d( fabs(y0_in[i]), fabs(y2_in[i]) ) ;
    diff = y2_in[i] - y2hat_in[i];
    err_it += (diff/sc)*(diff/sc);
  }
  return sqrt(err_it/( (double) dof));
}

int Cash_Karp::solve(double t_start, double t_end)
{
  int i, count;
  double h_new;
  hfull = h_init;
  t_it = t_start;
  dense_evals = 0;
  fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
  for ( i = 0; i < dof; i++)
  {
    y0[i] = y_init[i];
    fwrite(&(y0[i]), sizeof(double), 1, out_file_ptr);
  }

  eval(t_it, y0, k_new);

  count = 1;
  while (t_it < t_end)
  {
    count++;

    if ( (t_end-t_it) < h_new )
    {
      h_new = h_select(y0, t_end-t_it);
    }
    else
    {
      h_new = h_select(y0, hfull);
    }
    eval(t_it, yhat2, k_new);

    fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
    for ( i = 0; i < dof; i++)
    {
      y0[i] = yhat2[i];
      fwrite(&(y0[i]), sizeof(double), 1, out_file_ptr);
    }
    hfull = h_new;
  }
  aysml_gen(prefix, count, dof+1);
  return 1;
}

void Cash_Karp::dense_interp(double del_t)
{
  double a_0, a_1, pow_it, y_it, t_old, t_diff, theta_it, theta_acc;
  int signum;
  int n = 4;
  int i, j;
  t_old = t_it - 2.0*h;

  for ( i = 0; i < dof; i++) // gonna have to do the whole thing in the loop, otherwise, we're going to need to store multiple solves.
  {
    a_0 = y0[i];
    a_1 = 2*h*k_old[i];
    gsl_vector_set(b_vec, 0, yhat1[i] - a_0 - 0.5*a_1);
    gsl_vector_set(b_vec, 1, yhat2[i] - a_0 - a_1);
    gsl_vector_set(b_vec, 2, 2.0*h*k1[i] - a_1);
    gsl_vector_set(b_vec, 3, 2.0*h*k_new[i] - a_1);

    pow_it = 0.5;
    for ( j = 0; j < 4; j++)
    {
      pow_it *= pow_it;
      gsl_matrix_set(Theta, 0, j, pow_it);
      gsl_matrix_set(Theta, 1, j, 1.0);
      gsl_matrix_set(Theta, 2, j, ( (double) j + 2 )*pow_it);
      gsl_matrix_set(Theta, 3, j, ( (double) j + 2 ));
    }

    gsl_linalg_LU_decomp(Theta, perm_it, &signum);
    gsl_linalg_LU_solve(Theta, perm_it, b_vec, a_vec);
    gsl_matrix_set(a_full, i, 0, a_0);
    gsl_matrix_set(a_full, i, 1, a_1);
    gsl_matrix_set(a_full, i, 2, gsl_vector_get(a_vec, 0) );
    gsl_matrix_set(a_full, i, 3, gsl_vector_get(a_vec, 1) );
    gsl_matrix_set(a_full, i, 4, gsl_vector_get(a_vec, 2) );
    gsl_matrix_set(a_full, i, 5, gsl_vector_get(a_vec, 3) );
  }

  while (t_dense < t_it)
  {
    dense_evals++;
    t_diff = t_dense - t_old;
    theta_it = t_diff/(2.0*h);

    fwrite(&(t_dense), sizeof(double), 1, out_file_dense_ptr);
    for ( i = 0; i < dof; i++)
    {
      theta_acc = theta_it;
      y_it = gsl_matrix_get( a_full, i, 0);
      for ( j = 1; j < 6; j++)
      {
        y_it += theta_acc*gsl_matrix_get(a_full, i, j);
        theta_acc *= theta_acc;
      }
      fwrite(&(y_it), sizeof(double), 1, out_file_dense_ptr);
    }

    t_dense += del_t;
  }

}

int Cash_Karp::dense_solve(double t_start, double t_end, double del_t)
{
  int i, count;
  double h_new;
  hfull = h_init;
  t_it = t_start;
  t_dense = t_start;
  fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
  for ( i = 0; i < dof; i++)
  {
    y0[i] = y_init[i];
    fwrite(&(y0[i]), sizeof(double), 1, out_file_ptr);
  }

  eval(t_it, y0, k_new);

  count = 1;
  while (t_it < t_end)
  {
    count++;
    if ( (t_end-t_it) < h_new )
    {
      h_new = h_select(y0, t_end-t_it);
    }
    else
    {
      h_new = h_select(y0, hfull);
    }
    eval(t_it, yhat2, k_new);
    dense_interp(del_t);

    fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
    for ( i = 0; i < dof; i++)
    {
      y0[i] = yhat2[i];
      fwrite(&(y0[i]), sizeof(double), 1, out_file_ptr);
    }
    hfull = h_new;
  }
  aysml_gen(prefix, count, 3);
  aysml_gen(prefix_dense, dense_evals, 3);

  return 1;
}
