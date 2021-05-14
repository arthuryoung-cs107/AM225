#ifndef ADM_HH
#define ADM_HH

#include "AYmat.hh"

class ADM
{
  public:
    ADM(AYmat * X_);
    ~ADM();

    const int M, N;
    double mu, lambda;

    AYmat * X00;
    AYmat * L;
    AYmat * S;
    AYmat * X_k;
    AYmat * Sigma_work;

    AYmat * workMN;

    gsl_matrix * V_gsl;
    gsl_matrix * workNN_gsl;

    gsl_vector * workN_gsl;
    gsl_vector * s_gsl;

    void solve(bool verbose = false);
    void write_out(char prefix[]);

  private:
    void SVT();
    void shrink(AYmat * X_in, double tau, bool flag=false);
    double err_eval(AYmat * X, AYmat * L, AYmat * S);
    double max(double a, double b) { if (a > b) {return a;} else {return b;} }
    double min(double a, double b) { if (a < b) {return a;} else {return b;} }
};

#endif
