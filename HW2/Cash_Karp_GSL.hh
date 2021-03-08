#ifndef Cash_Karp_GSL_HH  /* Include guard */
#define Cash_Karp_GSL_HH

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


class Cash_Karp_GSL
{
  public:
    char prefix[100];
    char prefix_dense[150];
    char specfile[200];
    char specfile_dense[250];
    FILE * out_file_ptr;
    FILE * out_file_dense_ptr;

    const int dof, dims;
    int evals, dense_evals;
    const double facmax, facmin, fac;
    double a_tol, r_tol, h_init, hfull, h, t_it, t_dense;

    gsl_matrix *k1, *k2, *k3, *k4, *k5, *k6, *k_old, *k_new;
    gsl_matrix *w_it, *y_init, *y0, *y1, *y2, *yw, *yhat1, *yhat2;
    gsl_matrix * work1, *work2;

    double *** a_full;
    gsl_matrix * Theta;
    gsl_vector *b_vec, *a_vec;
    gsl_permutation * perm_it;

    Cash_Karp_GSL(int dof_in, int dims_in);
    ~Cash_Karp_GSL();
    double err_est(gsl_matrix * y0_in, gsl_matrix * y2_in, gsl_matrix * y2hat_in);

    virtual void eval(double time, gsl_matrix * y_in, gsl_matrix * y_out) = 0;
    virtual void write_out() = 0;

    int dense_solve(double t_start, double t_end, double del_t);

  private:
    void extrap(gsl_matrix * yw_in , gsl_matrix * y_in1, gsl_matrix * y_in2, gsl_matrix * y_out1, gsl_matrix * y_out2);
    void dense_interp(double del_t);
    double h_select(gsl_matrix * y_in, double h_in);
    void step(gsl_matrix * y_in, gsl_matrix * y_k, double t);
    inline double min_d(double a,double b) {return a<b?a:b;}
    inline double max_d(double a,double b) {return a>b?a:b;}
};

#endif
