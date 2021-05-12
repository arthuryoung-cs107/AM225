#include "ADM.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

ADM::ADM(AYmat * X_): M(X_->M), N(X_->N), X00(X_->copy_gen()),
L(new AYmat(M, N)), S(new AYmat(M, N)), X_k(new AYmat(M, N)),
Sigma_work(new AYmat(N, N)),
workMN(new AYmat(M, N)),
V_gsl(gsl_matrix_alloc(N, N)), s_gsl(gsl_vector_alloc(N)),
workNN_gsl(gsl_matrix_alloc(N, N)), workN_gsl(gsl_vector_alloc(N))
{
  Sigma_work->init_0();
  L->GSL_init();
  S->GSL_init();
}

ADM::~ADM()
{
  delete X00;
  delete L;
  delete S;
  delete X_k;
  delete Sigma_work;

  gsl_matrix_free(V_gsl);
  gsl_matrix_free(workNN_gsl);

  gsl_vector_free(s_gsl);
  gsl_vector_free(workN_gsl);
}

void ADM::solve(bool verbose)
{
  int count;
  int max_count = 1000;
  double tol;
  char prefix[200];

  mu = ((double) M)*((double) N)/(4.0*X00->norm_1());
  lambda = 1.0/(sqrt(max((double) M, (double) N)));
  tol = 1e-7*X00->norm_frob();

  count = 0;

  L->init_0();
  S->init_0();
  X_k->init_0();
  while ( (err_eval(X00, L, S) > tol ) && (count < max_count) )
  {
    SVT();
    X_k->copy_set(S);
    X00->add(L, S, -1.0, (1.0/mu));
    shrink(S, lambda/mu);

    S->copy_set(workMN);
    X00->add(L, workMN, -1., -1.);
    X_k->add(workMN, X_k, mu, 0.0 );
    count++;

    if (verbose)
    {
      // if (count%10 == 0)
      if ((count+1)%2 == 0)
      // if (count < 10)
      {
        memset(prefix, 0, 199);
        snprintf(prefix, 200, "./dat_dir/test2_it%d", count);
        write_out(prefix);
      }

      printf("iteration %d: error = %f \n", count, err_eval(X00, L, S));
    }

  }
}

void ADM::write_out(char prefix[])
{
  char nameS[200];
  memset(nameS, 0, 199);
  snprintf(nameS, 200, "%s_S_out", prefix);
  char nameL[200];
  memset(nameL, 0, 199);
  snprintf(nameL, 200, "%s_L_out", prefix);

  S->fprintf_mat(nameS);
  L->fprintf_mat(nameL);

}

void ADM::SVT()
{
  int i;
  X_k->copy_set(workMN);
  X00->add(S, workMN, -1.0, (1.0/mu));
  workMN->svd(s_gsl, V_gsl, workN_gsl);

  for ( i = 0; i < N; i++) Sigma_work->set(i, i, gsl_vector_get(s_gsl, i));

  shrink(Sigma_work, 1.0/mu, true);
  Sigma_work->GSL_init();

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, (double) 1.0, Sigma_work->A_gsl, V_gsl, (double) 0.0, workNN_gsl);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, (double) 1.0, workMN->A_gsl, workNN_gsl, (double) 0.0, L->A_gsl);
  L->GSL_send();
}

void ADM::shrink(AYmat * X_in, double tau, bool flag)
{
  double val;
  if (flag)
  {
    for ( int i = 0; i < N; i++)
    {
      val = max(abs(X_in->get(i, i))-tau, 0.0);
      if (X_in->get(i, i) < 0) val*=-1.0;
      X_in->set(i, i, val);
    }
  }
  else
  {
    for ( int i = 0; i < M; i++)
    {
      for ( int j = 0; j < N; j++)
      {
        val = max(abs(X_in->get(i, j))-tau, 0.0);
        if (X_in->get(i, j) < 0) val*=-1.0;
        X_in->set(i, j, val);
      }
    }
  }
}

double ADM::err_eval(AYmat * X_in, AYmat * L_in, AYmat * S_in)
{
    S_in->copy_set(workMN);
    X_in->add(L_in, workMN, -1., -1.);
    return workMN->norm_frob();
}
