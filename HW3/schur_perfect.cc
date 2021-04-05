#include <cstdio>
#include "schur_perfect.hh"
#include "blas.h"

/** Construct given the width of each subdomain.
 * \param[in] n_ : The width of each subdomain in units. */
schur_perfect::schur_perfect(square_specs * S) : conj_grad(S->n_glue),
n_sqrs(S->n_sqrs), n_glue(S->n_glue), nn_glue(S->nn_glue), N(S->N), h(1./((double) N+1)), ih2(1/(h*h)),
n_vec(S->n_vec), nn_vec(S->nn_vec), coords(S->coords),
mat_indices(S->mat_indices), vec_indices(S->vec_indices),
grids((poisson_fft **) malloc((n_sqrs-1)*(sizeof(poisson_fft*)))),
f_full(new double[N*N]), v_sol(new double[N*N]), v_sol2(new double[N*N]),
fk(new double[N*N]), A_glue(S->A_glue), A_glueT(S->A_glueT),
Minv_mat(new double[n_glue*n_glue]), A_check(new double[n_glue*n_glue]),
fk_full((double **) malloc((n_sqrs)*(sizeof(double*))))
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
  int i, j;
  for ( i = 0; i < n_sqrs-1; i++)
  {
    delete [] grids[i];
  }
}

void schur_perfect::init(const std::function<double(double,double)>& f)
{
  int i, j;
  double x_loc, y_loc;
  for(int j=0;j<N;j++)
  {
      y_loc=( (double) j+1)*h;
      for(int i=0;i<N;i++) {
          x_loc=( (double) i+1)*h;
          f_full[i+N*j] = f(x_loc, y_loc);
          fk[mat_indices[j][i]] = f(x_loc, y_loc);
      }
  }
}

/** Solve Poisson using the Schur complement method.
 * \param[in] f : A function handle to the source term. */
// void schur_perfect::solve_S(const std::function<double(double,double)>& f)
void schur_perfect::solve_S()
{
  int i, j, k, nn, ng, inc;
  double x_loc, y_loc, alpha, beta;
  char trans;

  inc = 1;
  ng = n_glue;

  for ( i = 0; i < n_glue; i++) b[i] = fk_full[n_sqrs-1][i];

  alpha = -1.0*ih2;
  beta = 1.0;
  trans = 't';

  for ( k = 0; k < n_sqrs-1; k++)
  {
    nn = nn_vec[k];
    for ( j = 0; j < nn_vec[k]; j++) grids[k]->f[j] = fk_full[k][j];

    grids[k]->solve();
    dgemv_(&trans, &nn, &ng, &alpha, A_glueT[k][0], &nn, grids[k]->v, &inc, &beta, b, &inc);
  }


  init_preconditioning();
  conj_grad::solve_pre(1);

  ////////////////////////////////
  // debugging zone
  ////////////////////////////////

  // alpha = 1.0;
  // beta = 0.0;
  // trans = 'n';
  //
  // char name[200];
  // memset(name, 0, 199);
  // snprintf(name, 100, "./dat_dir/prob5_Poisson_Schur_conjcheck");
  // for ( i = 0; i < N*N; i++) v_sol2[i] = v_sol[i] = 0;
  //
  // for ( j = 0; j < n_glue; j++)
  // {
  //   v_sol2[j+(N*N - n_glue)] = x[j];
  // }
  // for ( i = 0; i < N; i++)
  // {
  //   for ( j = 0; j < N; j++)
  //   {
  //     v_sol[i*N + j] = v_sol2[mat_indices[i][j]];
  //   }
  // }
  // printf("done\n");
  // output_solution(name, v_sol);
  // getchar();

  ////////////////////////////////
  // debugging zone
  ////////////////////////////////




  alpha = -1.0*ih2;
  beta = 1.0;
  trans = 'n';
  i = 0;
  for ( k = 0; k < n_sqrs-1; k++)
  {
    nn = nn_vec[k];
    for ( j = 0; j < nn; j++) grids[k]->f[j] = fk_full[k][j];
    dgemv_(&trans, &nn, &ng, &alpha, A_glueT[k][0], &nn, x, &inc, &beta, grids[k]->f, &inc);
    grids[k]->solve();

    for ( j = 0; j < nn; j++)
    {
      v_sol2[j+i] = grids[k]->v[j];
    }
    i += nn;
  }
  for ( j = 0; j < n_glue; j++)
  {
    v_sol2[j+i] = x[j];
  }
  for ( i = 0; i < N; i++)
  {
    for ( j = 0; j < N; j++)
    {
      v_sol[i*N + j] = v_sol2[mat_indices[i][j]];
    }
  }

}

/** Perform matrix-vector multiplication with the Schur complement matrix S.
 * \param[in] in  : The input vector.
 * \param[in] out : The output vector, after multiplying by S. */
void schur_perfect::mul_A(double* in, double* out)
{
  int k, nn, i;
  double alpha, beta;
  char trans = 'n';
  int inc = 1;
  int ng = n_glue;

  for ( i = 0; i < n_glue; i++) out[i] = 0.0;

  for ( k = 0; k < n_sqrs-1; k++)
  {
    nn = nn_vec[k];

    alpha = 1.0*ih2;
    beta = 0.0;
    trans = 'n';

    dgemv_(&trans, &nn, &ng, &alpha, A_glueT[k][0], &nn, in, &inc, &beta, grids[k]->f, &inc);
    grids[k]->solve();

    alpha = -1.0;
    beta = 1.0;
    trans = 't';
    dgemv_(&trans, &nn, &ng, &alpha, A_glueT[k][0], &nn, grids[k]->v, &inc, &beta, out, &inc);
  }
  alpha = 1.0;
  beta = 1.0;
  trans = 'n';
  // matrix vector multiplication with Aglue,glue
  dgemv_(&trans, &ng, &ng, &alpha, A_glueT[n_sqrs-1][0], &ng, in, &inc, &beta, out, &inc);

  for ( i = 0; i < n_glue; i++) out[i] *= ih2;
}

void schur_perfect::output_solution( char prefix[], double * v_in) {
    int n = N;
    const int ne=n+2;
    double *fld=new double[(n+2)*(n+2)],*f2=fld;

    // Set first row to zero
    while(f2<fld+ne) *(f2++)=0.;

    // Copy field contents into output array, padding the start and end entries
    // with zeros
    for(int j=0;j<n;j++) {
        *(f2++)=0;
        for(int i=0;i<n;i++) *(f2++)=v_in[i+j*n];
        *(f2++)=0;
    }

    // Set last row to zero, call output routine, and free output array
    while(f2<fld+ne*ne) *(f2++)=0.;
    write_out(prefix, ne, fld);
    delete [] fld;
}
void schur_perfect::write_out( char prefix[], int n, double * M)
{
  int i, j;
  FILE * out_file_ptr;
  char specfile[300];

  memset(specfile, 0, 299);
  snprintf(specfile, 200, "%s.aydat", prefix);
  out_file_ptr = fopen(specfile, "wb");

  for ( j = 0; j < n; j++)
  {
    for ( i = 0; i < n; i++) // walk
    {
      fwrite(&(M[i+n*j]), sizeof(double), 1, out_file_ptr);
    }
  }
  aysml_gen(prefix, n, n);

  fclose(out_file_ptr);
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
    for ( j = 0; j < n_glue; j++) A_check[i*n_glue + j] = 1.0/(work2[j]);
  }
  delete []work1;
  delete []work2;
}

void schur_perfect::M_inv(double *in,double *out)
{
  double alpha = 1.0;
  double beta = 0.0;
  char trans = 'n';
  int inc = 1;
  int ng = n_glue;
  dgemv_(&trans, &ng, &ng, &alpha, Minv_mat, &ng, in, &inc, &beta, out, &inc);
  // copy(in,out);

}
