#include <cstdio>
#include "schur_perfect.hh"
#include "blas.h"

/** Construct given the width of each subdomain.
 * \param[in] n_ : The width of each subdomain in units. */
schur_perfect::schur_perfect(square_specs * S) : conj_grad(S->n_glue),
n_sqrs(S->n_sqrs), n_glue(S->n_glue), nn_glue(S->nn_glue), N(S->N), h(1./((double) N+1)), ih2(1/(h*h)),
n_vec(S->n_vec), nn_vec(S->nn_vec),
mat_indices(S->mat_indices), vec_indices(S->vec_indices),
grids((poisson_fft **) malloc((n_sqrs-1)*(sizeof(poisson_fft*)))),
f_full(new double[N*N]), v_sol(new double[N*N]), Minv_mat(new double[n_glue*n_glue]),
fk_full((double **) malloc((n_sqrs)*(sizeof(double*)))),
fk(new double[N*N]), A_glue(S->A_glue)
{
  int i, j;
  for ( i = 0; i < n_sqrs-1; i++)
  {
    grids[i] = new poisson_fft(n_vec[i], h); // last one is solved by conjugate gradient
  }

  j = 0;
  for ( i = 0; i < n_sqrs; i++)
  {
    fk_full[i] = fk + j;
    j += nn_vec[i];
  }
}

schur_perfect::~schur_perfect()
{

}

/** Solve Poisson using the Schur complement method.
 * \param[in] f : A function handle to the source term. */
void schur_perfect::solve_S(const std::function<double(double,double)>& f)
{
  int i, j, k, n, ng, inc;
  double x_loc, y_loc, alpha, beta;
  char trans;

  inc = 1;
  ng = n_glue;

  for(int j=0;j<N;j++)
  {
      y_loc=( (double) j+1)*h;
      for(int i=0;i<N;i++) {
          x_loc=( (double) i+1)*h;
          f_full[i+N*j] = f(x_loc, y_loc);
          fk[mat_indices[j][i]] = f(x_loc, y_loc);
      }
  }

  for ( i = 0; i < n_glue; i++) b[i] = fk_full[n_sqrs-1][i];

  alpha = -1.0*ih2;
  beta = 1.0;
  trans = 't';

  for ( k = 0; k < n_sqrs-1; k++)
  {
    n = n_vec[k];
    for ( j = 0; j < nn_vec[k]; j++) grids[k]->f[j] = fk_full[k][j];

    grids[k]->solve();
    dgemv_(&trans, &n, &ng, &alpha, A_glue[k][0], &n, grids[k]->v, &inc, &beta, b, &inc);
  }

  init_preconditioning();
  conj_grad::solve_pre(1); // vglue solved, stored in x


  alpha = -1.0;
  beta = 1.0;
  trans = 'n';
  i = 0;
  for ( k = 0; k < n_sqrs-1; k++)
  {
    n = n_vec[k];
    for ( j = 0; j < nn_vec[k]; j++) grids[k]->f[j] = fk_full[k][j];
    dgemv_(&trans, &n, &ng, &alpha, A_glue[k][0], &n, x, &inc, &beta, grids[k]->f, &inc);
    grids[k]->solve();

    for ( j = 0; j < nn_vec[k]; j++)
    {
      v_sol[N*mat_indices[i+j][0] + mat_indices[i+j][1]] = grids[k]->v[j];
    }
    i += nn_vec[k];
  }
  for ( j = 0; j < nn_vec[k]; j++)
  {
    v_sol[N*mat_indices[i+j][0] + mat_indices[i+j][1]] = x[j];
  }
}

/** Perform matrix-vector multiplication with the Schur complement matrix S.
 * \param[in] in  : The input vector.
 * \param[in] out : The output vector, after multiplying by S. */
void schur_perfect::mul_A(double* in, double* out)
{
  int k, n, i;
  double alpha, beta;
  char trans = 'n';
  int inc = 1;
  int ng = n_glue;

  for ( i = 0; i < n_glue; i++) out[i] = 0.0;

  for ( k = 0; k < n_sqrs-1; k++)
  {
    n = n_vec[k];

    alpha = 1.0*ih2;
    beta = 0.0;
    trans = 'n';

    dgemv_(&trans, &n, &ng, &alpha, A_glue[k][0], &n, in, &inc, &beta, grids[k]->f, &inc);
    grids[k]->solve();

    alpha = -1.0;
    beta = 1.0;
    trans = 't';
    dgemv_(&trans, &n, &ng, &alpha, A_glue[k][0], &n, grids[k]->v, &inc, &beta, out, &inc);
  }

  alpha = 1.0;
  beta = 1.0;
  trans = 'n';
  // matrix vector multiplication with Aglue,glue
  dgemv_(&trans, &ng, &ng, &alpha, A_glue[n_sqrs-1][0], &ng, in, &inc, &beta, out, &inc);

  for ( i = 0; i < n_glue; i++) out[i] *= ih2;

}

/** Print the solution, padding with zeros for the Dirichlet boundary conditions. */
void schur_perfect::print_solution()
{

}

void schur_perfect::init_preconditioning()
{
  int i, j;
  double * work1 = new double[n_glue];
  double * work2 = new double[n_glue];

  for ( i = 0; i < n_glue*n_glue; i++) Minv_mat[i] = 0 ;
  for ( i = 0; i < n_glue; i++)
  {
    for ( j = 0; j < n_glue; j++) work1[j] = 0.0;
    work1[i] = 1.0;
    mul_A(work1, work2);
    Minv_mat[i*n_glue + i] = 1.0/(work2[i]);
  }
  delete []work1;
  delete []work2;
}

void schur_perfect::M_inv(double *in,double *out)
{
  // double alpha = 1.0;
  // double beta = 0.0;
  // char trans = 'n';
  // int inc = 1;
  // int ng = n_glue;
  // dgemv_(&trans, &ng, &ng, &alpha, Minv_mat, &ng, in, &inc, &beta, out, &inc);
  copy(in,out);

}
