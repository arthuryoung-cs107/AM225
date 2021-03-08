#include <cstdio>
#include <cmath>

#include "Cash_Karp_GSL.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

Cash_Karp_GSL::Cash_Karp_GSL(int dof_in, int dims_in): dof(dof_in), dims(dims_in), evals(0), h_init(0.01), fac(0.9), facmax(3.0), facmin(1.0/3.0),
k1(gsl_matrix_alloc(dof, dims)), k2(gsl_matrix_alloc(dof, dims)), k3(gsl_matrix_alloc(dof, dims)), k4(gsl_matrix_alloc(dof, dims)), k5(gsl_matrix_alloc(dof, dims)), k6(gsl_matrix_alloc(dof, dims)),
k_old(gsl_matrix_alloc(dof, dims)), k_new(gsl_matrix_alloc(dof, dims)),
w_it(gsl_matrix_alloc(dof, dims)), y_init(gsl_matrix_alloc(dof, dims)), y0(gsl_matrix_alloc(dof, dims)), y1(gsl_matrix_alloc(dof, dims)), y2(gsl_matrix_alloc(dof, dims)), yw(gsl_matrix_alloc(dof, dims)), yhat1(gsl_matrix_alloc(dof, dims)), yhat2(gsl_matrix_alloc(dof, dims)),
work1(gsl_matrix_alloc(dof, dims)), work2(gsl_matrix_alloc(dof, dims)),
a_full(d3tensor(0, dof-1, 0, dims-1, 0, 5)),
Theta(gsl_matrix_alloc(4, 4)), b_vec(gsl_vector_alloc(4)), a_vec(gsl_vector_alloc(4)), perm_it(gsl_permutation_alloc(4))
{}

Cash_Karp_GSL::~Cash_Karp_GSL()
{
  gsl_matrix_free(k1);
  gsl_matrix_free(k2);
  gsl_matrix_free(k3);
  gsl_matrix_free(k4);
  gsl_matrix_free(k5);
  gsl_matrix_free(k6);
  gsl_matrix_free(k_old);
  gsl_matrix_free(k_new);
  gsl_matrix_free(w_it);
  gsl_matrix_free(y_init);
  gsl_matrix_free(y0);
  gsl_matrix_free(y1);
  gsl_matrix_free(y2);
  gsl_matrix_free(yw);
  gsl_matrix_free(yhat1);
  gsl_matrix_free(yhat2);

  gsl_matrix_free(work1);
  gsl_matrix_free(work2);

  free_d3tensor(a_full, 0, dof-1, 0, dims-1, 0, 5);
  gsl_matrix_free(Theta);
  gsl_vector_free(b_vec);
  gsl_vector_free(a_vec);
  gsl_permutation_free(perm_it);
}

void Cash_Karp_GSL::step(gsl_matrix * y_in, gsl_matrix * y_k, double t)
{
  int i, out;
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

  out = gsl_matrix_memcpy(w_it, y_in);
  AY_GSLmatrix_add(w_it, k1, work1, h*a21);
  eval(t + c2*h, w_it, k2);

  out = gsl_matrix_memcpy(w_it, y_in);
  AY_GSLmatrix_add(w_it, k1, work1, h*a31);
  AY_GSLmatrix_add(w_it, k2, work1, h*a32);
  eval(t + c3*h, w_it, k3);

  out = gsl_matrix_memcpy(w_it, y_in);
  AY_GSLmatrix_add(w_it, k1, work1, h*a41);
  AY_GSLmatrix_add(w_it, k2, work1, h*a42);
  AY_GSLmatrix_add(w_it, k3, work1, h*a43);
  eval(t + c4*h, w_it, k4);

  out = gsl_matrix_memcpy(w_it, y_in);
  AY_GSLmatrix_add(w_it, k1, work1, h*a51);
  AY_GSLmatrix_add(w_it, k2, work1, h*a52);
  AY_GSLmatrix_add(w_it, k3, work1, h*a53);
  AY_GSLmatrix_add(w_it, k4, work1, h*a54);
  eval(t + c5*h, w_it, k5);

  out = gsl_matrix_memcpy(w_it, y_in);
  AY_GSLmatrix_add(w_it, k1, work1, h*a61);
  AY_GSLmatrix_add(w_it, k2, work1, h*a62);
  AY_GSLmatrix_add(w_it, k3, work1, h*a63);
  AY_GSLmatrix_add(w_it, k4, work1, h*a64);
  AY_GSLmatrix_add(w_it, k5, work1, h*a65);
  eval(t + c6*h, w_it, k6);

  // add em up
  out = gsl_matrix_memcpy(y_k, y_in);
  AY_GSLmatrix_add(y_k, k1, work1, h*b1);
  AY_GSLmatrix_add(y_k, k2, work1, h*b2);
  AY_GSLmatrix_add(y_k, k3, work1, h*b3);
  AY_GSLmatrix_add(y_k, k4, work1, h*b4);
  AY_GSLmatrix_add(y_k, k5, work1, h*b5);
  AY_GSLmatrix_add(y_k, k6, work1, h*b6);
}

void Cash_Karp_GSL::extrap(gsl_matrix * yw_in , gsl_matrix * y_in1, gsl_matrix * y_in2, gsl_matrix * y_out1, gsl_matrix * y_out2)
{
  gsl_matrix_memcpy(y_out1, y_in1);
  gsl_matrix_memcpy(y_out2, y_in2);
  gsl_matrix_memcpy(work1, y_in2);
  gsl_matrix_sub(work1, yw_in);
  AY_GSLmatrix_add(y_out1, work1, work2, (1.0)/(62.0));
  AY_GSLmatrix_add(y_out2, work1, work2, (1.0)/(31.0));
}

double Cash_Karp_GSL::h_select(gsl_matrix * y_in, double h_in)
{
  int i;
  double err, h_new, h_it;
  h_it = h_in;
  gsl_matrix_memcpy(k_old, k_new);
  do
  {
    gsl_matrix_memcpy(k1, k_old);
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

double Cash_Karp_GSL::err_est(gsl_matrix * y0_in, gsl_matrix * y2_in, gsl_matrix * y2hat_in)
{
  double sc, diff;
  double err_it = 0;
  gsl_matrix_memcpy(work1, y2_in);
  gsl_matrix_sub(work1, y2hat_in);
  for (int i = 0; i < dof; i++)
  {
    for (int j = 0; j < dims; j++)
    {
      sc = a_tol + r_tol*max_d( fabs(gsl_matrix_get(y0_in, i, j) ), fabs(gsl_matrix_get(y2_in, i, j)) ) ;
      err_it += (gsl_matrix_get(work1, i, j)/sc)*(gsl_matrix_get(work1, i, j)/sc);
    }
  }
  return sqrt(err_it/( (double) dof));
}

void Cash_Karp_GSL::dense_interp(double del_t)
{
  double a_0, a_1, pow_it, y_it, t_old, t_diff, theta_it, theta_acc;
  int signum;
  int n = 4;
  int i, j, k;
  t_old = t_it - 2.0*h;

  for ( i = 0; i < dof; i++) // gonna have to do the whole thing in the loop, otherwise, we're going to need to store multiple solves.
  {
    for ( j = 0; j < dims; j++)
    {
      a_0 = gsl_matrix_get(y0, i, j);
      a_1 = 2*h*gsl_matrix_get(k_old, i, j);
      gsl_vector_set(b_vec, 0, gsl_matrix_get(yhat1, i, j) - a_0 - 0.5*a_1);
      gsl_vector_set(b_vec, 1, gsl_matrix_get(yhat2, i, j) - a_0 - a_1);
      gsl_vector_set(b_vec, 2, 2.0*h*gsl_matrix_get(k1, i, j) - a_1);
      gsl_vector_set(b_vec, 3, 2.0*h*gsl_matrix_get(k_new, i, j) - a_1);

      pow_it = 0.5;
      for ( k = 0; k < 4; k++)
      {
        pow_it *= pow_it;
        gsl_matrix_set(Theta, 0, k, pow_it);
        gsl_matrix_set(Theta, 1, k, 1.0);
        gsl_matrix_set(Theta, 2, k, ( (double) k + 2 )*pow_it);
        gsl_matrix_set(Theta, 3, k, ( (double) k + 2 ));
      }

      gsl_linalg_LU_decomp(Theta, perm_it, &signum);
      gsl_linalg_LU_solve(Theta, perm_it, b_vec, a_vec);

      a_full[i][j][0] = a_0;
      a_full[i][j][1] = a_1;
      a_full[i][j][2] = gsl_vector_get(a_vec, 0);
      a_full[i][j][3] = gsl_vector_get(a_vec, 1);
      a_full[i][j][4] = gsl_vector_get(a_vec, 2);
      a_full[i][j][5] = gsl_vector_get(a_vec, 3);
    }
  }

  while (t_dense < t_it)
  {
    dense_evals++;
    t_diff = t_dense - t_old;
    theta_it = t_diff/(2.0*h);

    fwrite(&(t_dense), sizeof(double), 1, out_file_dense_ptr);
    for ( i = 0; i < dof; i++)
    {
      for ( j = 0; j < dims; j++)
      {
        theta_acc = theta_it;
        y_it = a_full[i][j][0];
        for ( k = 1; k < 6; k++)
        {
          y_it += theta_acc*a_full[i][j][k];
          theta_acc *= theta_acc;
        }
        fwrite(&(y_it), sizeof(double), 1, out_file_dense_ptr);
      }
    }
    t_dense += del_t;
  }
}

int Cash_Karp_GSL::dense_solve(double t_start, double t_end, double del_t)
{
  out_file_ptr = fopen(specfile, "wb");
  out_file_dense_ptr = fopen(specfile_dense, "wb");
  int i, count;
  double h_new;
  hfull = h_init;
  t_it = t_start;
  t_dense = t_start;
  dense_evals = 0;
  gsl_matrix_memcpy(y0, y_init);
  write_out();
  eval(t_it, y0, k_new);
  count = 1;
  while (t_it < t_end)
  {
    count++;

    if (count%100 == 0)
    {
      printf("count: %d, time: %f out of %f\n", count, t_it, t_end);
    }

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

    gsl_matrix_memcpy(y0, yhat2);
    write_out();
    hfull = h_new;
  }
  fclose(out_file_dense_ptr);
  fclose(out_file_ptr);
  aysml_gen(prefix, count, dof*dims);
  aysml_gen(prefix_dense, dense_evals, dof*dims);

  return 1;
}
