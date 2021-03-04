#ifndef Cash_Karp_HH  /* Include guard */
#define Cash_Karp_HH

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


class Cash_Karp
{
  public:
    char prefix[100];
    char prefix_dense[150];
    char specfile[200];
    char specfile_dense[250];
    FILE * out_file_ptr;
    FILE * out_file_dense_ptr;

    const int dof;
    int evals, dense_evals;
    const double facmax, facmin, fac;
    double a_tol, r_tol, h_init, hfull, h, t_it, t_dense;
    double *k1, *k2, *k3, *k4, *k5, *k6, *k_old, *k_new;
    double *w_it, *y_init, *y0, *y1, *y2, *yw, *yhat1, *yhat2;

    gsl_matrix * Theta, *a_full;
    gsl_vector *b_vec, *a_vec;
    gsl_permutation * perm_it;


    Cash_Karp(int dof_in);
    ~Cash_Karp();
    double err_est(double * y0_in, double * y2_in, double * y2hat_in);

    virtual void eval(double time, double * y_in, double * y_out) = 0;

    int solve(double t_start, double t_end);
    int dense_solve(double t_start, double t_end, double del_t);

  private:
    void extrap(double * yw_in , double * y_in1, double * y_in2, double * y_out1, double * y_out2);
    void dense_interp(double del_t);
    double h_select(double * y_in, double h_in);
    void step(double * y_in, double * y_k, double t);
    inline double min_d(double a,double b) {return a<b?a:b;}
    inline double max_d(double a,double b) {return a>b?a:b;}
};

#endif
