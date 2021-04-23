#include "Ritz_Galerk_sphere.hh"
#include "blas.h"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

Ritz_Galerk_sphere::Ritz_Galerk_sphere(int n_) : conj_grad(((n_-1)*(n_-1))),
dof((n_-1)*(n_-1)),  n_full(n_), n(n_-1), n_q(10), xy(new double[2]), vw(new double[2]),
vw_coords(dmatrix(0, dof-1, 0, 1)), xy_coords(dmatrix(0, dof-1, 0, 1)),
h(2.0/((double) n_)),
Jac(new AYmat(2, 2)), a_vals(new AYmat(dof, 9)), Jac_inv(new AYmat(2, 2)),
dphi_i(new AYmat(2, 1)), dphi_j(new AYmat(2, 1)), work1(new AYmat(2, 1)), work2(new AYmat(2, 1)),
a_count(new int[dof]), a_indices(imatrix(0, dof-1, 0, 8)), q(new quadrat(n_q))
{
  int i, j;
  a_vals->init_0();

  for ( i = 0; i < dof; i++) a_count[i] = 0 ;

  for ( i = 0; i < n; i++)
  {
    for ( j = 0; j < n; j++)
    {
      xy_coords[n*(i) + j][0] = -1.0 + h*((double) (i + 1) );
      xy_coords[n*(i) + j][1] = -1.0 + h*((double) (j + 1) );
      map(xy_coords[n*(i) + j], vw_coords[n*(i) + j]);
    }
  }
  assemble_a();

}

Ritz_Galerk_sphere::~Ritz_Galerk_sphere()
{
  delete [] xy;
  delete [] vw;
  delete [] a_count;
  delete Jac;
  delete a_vals;
  delete Jac_inv;
  delete dphi_i;
  delete dphi_j;
  delete work1;
  delete work2;
  delete q;
  free_imatrix(a_indices, 0, dof-1, 0, 8);
  free_dmatrix(xy_coords, 0, dof-1, 0, 1);
  free_dmatrix(vw_coords, 0, dof-1, 0, 1);
}

/** Performs multiplication on a vector by the stiffness matrix. */
void Ritz_Galerk_sphere::mul_A(double *in,double *out)
{
  int i, j;
  double sum;
  for ( i = 0; i < dof; i++)
  {
    sum = 0;
    for ( j = 0; j < a_count[i]; j++)
    {
      sum += a_vals->get( i, j )*in[a_indices[i][j]];
    }
    out[i] = sum;
  }
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

void Ritz_Galerk_sphere::assemble_b(const std::function<double(double,double)>& f_source)
{
  int i, j, k;
  double acc, acc_r, x_c, y_c;

  for ( i = 0; i < dof; i++)
  {
    acc = 0;
    x_c = xy_coords[i][0];
    y_c = xy_coords[i][1];
    for ( j = 0; j < n_q; j++)
    {
      acc_r = 0;
      xy[1] = y_c + q->x[j]*h;
      for ( k = 0; k < n_q; k++)
      {
        xy[0] = x_c + q->x[k]*h;
        map(xy, vw);
        Jac_eval(xy, Jac);
        acc_r += q->w[k]*f_source(vw[0], vw[1])*phi_eval(q->x[k], q->x[j])*Jac_det(Jac);
      }
      acc+=q->w[j]*acc_r;
    }
    b[i] = acc*h*h;
  }
}

void Ritz_Galerk_sphere::assemble_a()
{
  double prod;
  int i, j, row, col;


  for ( i = 0; i < dof; i++) // identify all non-zero interactions.
  {
    row = i/n;
    col = i%n;

    a_vals->set(i, a_count[i], a_prod(i, i) );
    a_indices[i][a_count[i]] = i;
    a_count[i]++;

    if (row > 0) // if NOT on the TOP row, then there must be a node ABOVE
    {
      j = i - n; // index of above node
      a_vals->set(i, a_count[i], a_prod(i, j) );
      a_indices[i][a_count[i]] = j; // record an interaction
      a_count[i]++;

      if (col == n-1) // right edge?
      {
        j = i - n - 1; // index of NW node
        a_vals->set(i, a_count[i], a_prod(i, j) );
        a_indices[i][a_count[i]] = j; // record an interaction
        a_count[i]++;
      }
      else
      {
        if (col == 0) // left edge?
        {
          j = i - n + 1; // index of NE node
          a_vals->set(i, a_count[i], a_prod(i, j) );
          a_indices[i][a_count[i]] = j; // record an interaction
          a_count[i]++;
        }
        else // interior point.
        {
          j = i - n - 1; // index of NW node
          a_vals->set(i, a_count[i], a_prod(i, j) );
          a_indices[i][a_count[i]] = j; // record an interaction
          a_count[i]++;

          j = i - n + 1; // index of NE node
          a_vals->set(i, a_count[i], a_prod(i, j) );
          a_indices[i][a_count[i]] = j; // record an interaction
          a_count[i]++;
        }
      }
    }
    if (row < n-1) // if NOT on the BOTTOM row, then there must be a node BELOW
    {
      j = i + n; // index of below node
      a_vals->set(i, a_count[i], a_prod(i, j) );
      a_indices[i][a_count[i]] = j; // record an interaction
      a_count[i]++;

      if (col == n-1) // right edge?
      {
        j = i + n - 1; // index of SW node
        a_vals->set(i, a_count[i], a_prod(i, j) );
        a_indices[i][a_count[i]] = j; // record an interaction
        a_count[i]++;
      }
      else
      {
        if (col == 0) // left edge?
        {
          j = i + n + 1; // index of SE node
          a_vals->set(i, a_count[i], a_prod(i, j) );
          a_indices[i][a_count[i]] = j; // record an interaction
          a_count[i]++;
        }
        else // interior point.
        {
          j = i + n - 1; // index of SW node
          a_vals->set(i, a_count[i], a_prod(i, j) );
          a_indices[i][a_count[i]] = j; // record an interaction
          a_count[i]++;

          j = i + n + 1; // index of SE node
          a_vals->set(i, a_count[i], a_prod(i, j) );
          a_indices[i][a_count[i]] = j; // record an interaction
          a_count[i]++;
        }
      }
    }
    if (col > 0) // if NOT on left column, then there must be a node WEST
    {
      j = i - 1; // index of west node
      a_vals->set(i, a_count[i], a_prod(i, j) );
      a_indices[i][a_count[i]] = j; // record an interaction
      a_count[i]++;
    }
    if (col < n-1) // if NOT on right column, then there must be a node EAST
    {
      j = i + 1; // index of east node
      a_vals->set(i, a_count[i], a_prod(i, j) );
      a_indices[i][a_count[i]] = j; // record an interaction
      a_count[i]++;
    }
  }

}

double Ritz_Galerk_sphere::a_prod(int i_in, int j_in)
{
  int i, j, k;
  double acc_full, acc, acc_r, x_c, y_c, x_j, y_j, del_x, del_y;

  acc = 0;
  x_c = xy_coords[i_in][0];
  y_c = xy_coords[i_in][1];

  x_j = xy_coords[j_in][0];
  y_j = xy_coords[j_in][1];

  del_x = (x_c-x_j)/(h);
  del_y = (y_c-y_j)/(h);

  for ( j = 0; j < n_q; j++)
  {
    acc_r = 0;
    xy[1] = y_c + (q->x[j]*h);
    for ( k = 0; k < n_q; k++)
    {
      xy[0] = x_c + (q->x[k]*h);
      if ( (abs(xy[0] - x_j) < h) && (abs(xy[1] - y_j) < h) ) // overlap with phi_j
      {
        Jac_eval(xy, Jac);
        Jac_invert(Jac, Jac_inv);
        grad_phi_eval(q->x[k], q->x[j], dphi_i);
        grad_phi_eval(del_x+q->x[k], del_y+q->x[j], dphi_j);
        work1->mult_set(Jac_inv, dphi_i, 1.0, 0.0);
        work2->mult_set(Jac_inv, dphi_j, 1.0, 0.0);
        acc_r += q->w[k]*(work1->inner(work2))*Jac_det(Jac);
      }
    }
    acc+=q->w[j]*acc_r;
  }
  acc_full = acc*h*h;
  return  acc_full;
}

void Ritz_Galerk_sphere::Jac_invert(AYmat * Jac_in, AYmat * Jac_inverted)
{
  double det = Jac_det(Jac_in);
  Jac_inverted->set( 0, 0, Jac_in->get(1, 1)/det );
  Jac_inverted->set( 0, 1, -Jac_in->get(1, 0)/det );
  Jac_inverted->set( 1, 0, -Jac_in->get(0, 1)/det );
  Jac_inverted->set( 1, 1, Jac_in->get(0, 0)/det );
}

void Ritz_Galerk_sphere::grad_phi_eval(double x_loc, double y_loc, AYmat * grad_out)
{
  if (x_loc > 0)
  {
    if (y_loc > 0) // I
    {
      grad_out->set(0, 0, (-1.0 + y_loc)*(1.0/h));
      grad_out->set(1, 0, (-1.0 + x_loc)*(1.0/h));
    }
    else // IV
    {
      grad_out->set(0, 0, (-1.0 - y_loc)*(1.0/h));
      grad_out->set(1, 0, (1.0 - x_loc)*(1.0/h));
    }
  }
  else
  {
    if (y_loc > 0) // II
    {
      grad_out->set(0, 0, (1.0 - y_loc)*(1.0/h));
      grad_out->set(1, 0, (-1.0 - x_loc)*(1.0/h));
    }
    else  // III
    {
      grad_out->set(0, 0, (1.0 + y_loc)*(1.0/h));
      grad_out->set(1, 0, (1.0 + x_loc)*(1.0/h));
    }
  }
}

void Ritz_Galerk_sphere::write_out(char prefix_x[], char prefix_y[], char prefix_v[], char prefix_w[], char prefix_u[], int N)
{
    int i, j, k ;
    double acc;

    double * test_coords = new double[N];
    double ** xyvw = dmatrix(0, N*N-1, 0, 3);
    double *** sol_out = d3tensor(0, 4, 0, N-1, 0, N-1);

    double del = 2.0/( (double) N - 1 );

    for ( i = 0; i < N; i++) test_coords[i] = -1.0 + del*( (double) i);

    for ( i = 0; i < N; i++)
    {
      for ( j = 0; j < N; j++)
      {
        xyvw[i*N+j][0] = test_coords[i]; xyvw[i*N+j][1] = test_coords[j];
        map(xyvw[i*N+j], xyvw[i*N+j]+2);
        for ( k = 0; k < 4; k++) sol_out[k][i][j] = xyvw[i*N+j][k];

        if (i == 0 || i == N-1 || j == 0 || j == N-1) // if edge/corner node
        {
          sol_out[4][i][j] = 0.0;
        }
        else
        {
          acc = 0;
          for ( k = 0; k < dof; k++)
          {
            if ( abs(xyvw[i*N+j][0] - xy_coords[k][0]) < h && abs(xyvw[i*N+j][1] - xy_coords[k][1]) < h )
            {
              acc += x[k]*phi_eval( (xyvw[i*N+j][0] - xy_coords[k][0])/h , (xyvw[i*N+j][1] - xy_coords[k][1])/h );
            }
          }
          sol_out[4][i][j] = acc;
        }
      }
    }
    fprintf_matrix(sol_out[0], N, N, prefix_x);
    fprintf_matrix(sol_out[1], N, N, prefix_y);
    fprintf_matrix(sol_out[2], N, N, prefix_v);
    fprintf_matrix(sol_out[3], N, N, prefix_w);
    fprintf_matrix(sol_out[4], N, N, prefix_u);

    delete [] test_coords;
    free_dmatrix(xyvw, 0, N*N-1, 0, 3);
    free_d3tensor(sol_out, 0, 4, 0, N-1, 0, N-1);
}
