#ifndef Cash_Karp_HH  /* Include guard */
#define Cash_Karp_HH

class Cash_Karp
{
  public:
    char prefix[100];
    char specfile[200];
    FILE * out_file_ptr;

    const int dof;
    const double facmax, facmin, fac;
    double a_tol, r_tol, h_init, hfull, h, t_it;
    double *k1, *k2, *k3, *k4, *k5, *k6;
    double *w_it, *y_init, *y0, *y1, *y2, *yw, *yhat1, *yhat2;

    Cash_Karp(int dof_in);
    ~Cash_Karp();
    double err_est(double * y0_in, double * y2_in, double * y2hat_in);

    virtual void eval(double time, double * y_in, double * y_out) = 0;

    int solve(double t_start, double t_end);

  private:
    void extrap(double * yw_in , double * y_in1, double * y_in2, double * y_out1, double * y_out2);
    double h_select(double * y_in, double h_in);
    void step(double * y_in, double * y_k, double t);
    inline double min_d(double a,double b) {return a<b?a:b;}
    inline double max_d(double a,double b) {return a>b?a:b;}
};

#endif
