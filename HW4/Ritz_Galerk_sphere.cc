#include "Ritz_Galerk_sphere.hh"
#include "quadrat.hh"
#include "blas.h"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

// Ritz_Galerk_sphere::Ritz_Galerk_sphere(int n_) : dof((n_-1)*(n_-1)), conj_grad(dof), n_full(n_), n(n_-1), f(new double[dof]),  vw_coords(dmatrix(0, dof-1, 0, 1)), xy_coords(dmatrix(0, dof-1, 0, 1)), h(2.0/n_full)
Ritz_Galerk_sphere::Ritz_Galerk_sphere(int n_) : conj_grad(((n_-1)*(n_-1))),
dof((n_-1)*(n_-1)),  n_full(n_), n(n_-1),
f(new double[dof]), xy(new double[2]), vw(new double[2]),
vw_coords(dmatrix(0, dof-1, 0, 1)), xy_coords(dmatrix(0, dof-1, 0, 1)),
h(2.0/n_full), Jac(new AYmat(2, 2))
{
  int i, j;
  AYmat M1 = AYmat(n, n-1);


  for ( i = 0; i < n; i++)
  {
    for ( j = 0; j < n; j++)
    {
      xy_coords[n*(i) + j][0] = -1.0 + h*((double) (i + 1) );
      xy_coords[n*(i) + j][1] = -1.0 + h*((double) (j + 1) );
      map(xy_coords[n*(i) + j], vw_coords[n*(i) + j]);
    }
  }
}

Ritz_Galerk_sphere::~Ritz_Galerk_sphere()
{
  delete [] f;
  delete [] xy;
  delete [] vw;
  free_dmatrix(xy_coords, 0, dof-1, 0, 1);
  free_dmatrix(vw_coords, 0, dof-1, 0, 1);
}

/** Performs multiplication on a vector by the stiffness matrix. */
void Ritz_Galerk_sphere::mul_A(double *in,double *out)
{

}

void Ritz_Galerk_sphere::map(double * x_in, double * v_out)
{
  double x = x_in[0]; double y = x_in[1];
  v_out[0] = x*sqrt(1.0 - (y*y*0.5) );
  v_out[1] = y*sqrt(1.0 - (x*x*0.5) );
}

void Ritz_Galerk_sphere::Jac_eval(double * x_in, AYmat * mat_out)
{
  double x = x_in[0]; double y = x_in[1];
  double J11, J12, J21, J22;

  J11 = sqrt( 1.0 - (y*y*0.5));
  J12 = 0.5*x*(1.0/sqrt(1.0 - y*y*0.5))*(-y);
  J21 = 0.5*y*(1.0/sqrt(1.0 - x*x*0.5))*(-x);
  J22 = sqrt( 1.0 - (x*x*0.5));

  mat_out->set(0, 0, J11);
  mat_out->set(0, 1, J12);
  mat_out->set(1, 0, J21);
  mat_out->set(1, 1, J22);
}

double Ritz_Galerk_sphere::Jac_det(AYmat * Jac_in) {return ( (Jac_in->get(0, 0) * Jac_in->get(1, 1) ) - (Jac_in->get(0, 1) * Jac_in->get(1, 0)) );}

double Ritz_Galerk_sphere::phi_eval(double x_loc, double y_loc)
{
  if (x_loc > 0)
  {
    if (y_loc > 0) { return (1.0 - x_loc - y_loc + x_loc*y_loc); } // I
    else { return(1.0 - x_loc + y_loc - x_loc*y_loc); } // IV
  }
  else
  {
    if (y_loc > 0) { return (1.0 + x_loc - y_loc - x_loc*y_loc); } // II
    else { return (1.0 + x_loc + y_loc + x_loc*y_loc); }  // III
  }
}

/** Computes the source vector in the linear system, which is based
 * on multiplying the source function by the mass matrix. */
void Ritz_Galerk_sphere::assemble_b(const std::function<double(double,double)>& f_source)
{
  int i, j, k;
  double acc, acc_r, x_c, y_c;
  int n_q = 5;
  quadrat q(n_q);

  for ( i = 0; i < dof; i++)
  {
    acc = 0;
    x_c = xy_coords[i][0];
    y_c = xy_coords[i][1];
    for ( j = 0; j < n_q; j++)
    {
      acc_r = 0;
      xy[1] = y_c + q.x[j]*h;
      for ( k = 0; k < n_q; k++)
      {
        xy[0] = x_c + q.x[k]*h;
        map(xy, vw);
        Jac_eval(xy, Jac);
        acc_r += q.w[k]*f_source(vw[0], vw[1])*phi_eval(q.x[k], q.x[j])*Jac_det(Jac);
      }
      acc+=q.w[j]*acc_r;
    }
    b[i] = acc*h*h;
  }


}
