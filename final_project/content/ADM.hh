#ifndef ADM_HH
#define ADM_HH

#include "AYmat.hh"

class ADM
{
  public:
    const int M, N;
    double mu, lambda;

    AYmat * X00;
    AYmat * L;
    AYmat * S;
    AYmat * X_k;
    AYmat * Sigma_work;

    AYmat * workMN;

    gsl_matrix * V_gsl;
    gsl_vector * s_gsl;

    gsl_matrix * workNN_gsl;
    gsl_vector * workN_gsl;

    ADM(AYmat * X_);
    ~ADM();

    void solve();

  private:
    void SVT();
    void shrink(AYmat * X_in, double tau, int flag=1);
    double err_eval(AYmat * X, AYmat * L, AYmat * S);
    double max(double a, double b) { if (a > b) {return a;} else {return b;} }
    double min(double a, double b) { if (a < b) {return a;} else {return b;} }
};

#endif
