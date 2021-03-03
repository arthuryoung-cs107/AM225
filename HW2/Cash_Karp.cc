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
w_it(new double[dof]), y_init(new double[dof]), y0(new double[dof]), y1(new double[dof]), y2(new double[dof]), yw(new double[dof]), yhat1(new double[dof]), yhat2(new double[dof]) {}

Cash_Karp::~Cash_Karp()
{
  delete [] k1;
  delete [] k2;
  delete [] k3;
  delete [] k4;
  delete [] k5;
  delete [] k6;
  delete [] w_it;
  delete [] y_init;
  delete [] y0;
  delete [] y1;
  delete [] y2;
  delete [] yw;
  delete [] yhat1;
  delete [] yhat2;
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
  do
  {
    h = h_it;
    eval(t_it, y_in, k1);
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
  fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
  for ( i = 0; i < dof; i++)
  {
    y0[i] = y_init[i];
    fwrite(&(y0[i]), sizeof(double), 1, out_file_ptr);
  }


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

    fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
    for ( i = 0; i < dof; i++)
    {
      y0[i] = yhat2[i];
      fwrite(&(y0[i]), sizeof(double), 1, out_file_ptr);
    }
    hfull = h_new;
  }
  fclose(out_file_ptr);
  aysml_gen(prefix, count, 3);
  return 1;
}
