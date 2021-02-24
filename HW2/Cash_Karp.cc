#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <string.h>

#include "Cash_Karp.hh"

Cash_Karp::Cash_Karp(int dof_in): dof(dof_in),
atol(1e-14), rtol(1e-14),
w(new double[dof]), w_it(new double[dof]), k1(new double[dof]), k2(new double[dof]), k3(new double[dof]), k4(new double[dof]), k5(new double[dof]), k6(new double[dof]) {}

Cash_Karp::~Cash_Karp()
{
  delete [] w;
  delete [] w_it;
  delete [] k1;
  delete [] k2;
  delete [] k3;
  delete [] k4;
  delete [] k5;
  delete [] k6;
}

double



void Cash_Karp::step(double * y_k)
{
  int i;
  const double c2 = 0.2,
               c3 = 0.3,
               c4 = 0.6,
               c5 = 1.0,
               c6 = 0.875,

               b1 = 97.8835978835979e-003,
               b2 = 0.0,
               b3 = 402.576489533011e-003,
               b4 = 210.437710437710e-003,
               b5 = 0.0,
               b6 = 289.102202145680e-003,

               a21 = 0.2,
               a31 = 0.075,
               a32 = 0.225,
               a41 = 0.3,
               a42 = -0.9,
               a43 = 1.2,
               a51 = -203.703703703704e-003,
               a52 = 2.5,
               a53 = -2.59259259259259e+000,
               a54 = 1.29629629629630e+000,
               a61 = 29.4958043981481e-003,
               a62 = 341.796875000000e-003,
               a63 = 41.5943287037037e-003,
               a64 = 400.345413773148e-003,
               a65 = 61.7675781250000e-003;

  // assume k1 initialized to current function value
  for (i = 0; i < n; i++) w_it[i] = w[i] + h*a21*k1[i];
  eval(t + c2*h, w_it, k2);

  for (i = 0; i < n; i++) w_it[i] = w[i] + h*(a31*k1[i] + a32*k2[i]);
  eval(t + c3*h, w_it, k3);

  for (i = 0; i < n; i++) w_it[i] = w[i] + h*(a41*k1[i] + a42*k2[i] + a43*k3[i]);
  eval(t + c4*h, w_it, k4);

  for (i = 0; i < n; i++) w_it[i] = w[i] + h*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
  eval(t + c5*h, w_it, k5);

  for (i = 0; i < n; i++) w_it[i] = w[i] + h*(a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
  eval(t + c6*h, w_it, k6);

  // add em up
  for (i = 0; i < n; i++) y_k[i] = w[i] + h*( b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i]);
}

int Cash_Karp::solve(double t_start, double t_end, double ** outputs)
{

}
