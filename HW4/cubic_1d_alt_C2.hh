#ifndef CUBIC_1D_ALT_C2_HH
#define CUBIC_1D_ALT_C2_HH

#include "conj_grad.hh"
#include "quadrat.hh"

#include <cstdio>
#include <cmath>

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

/** Class for solving an elliptic PDE problem over the domain [1,2] using
 * piecewise cubic basis functions. */
class cubic_1d_alt_C2 : public conj_grad {
    public:
        /** The number of intervals to divide the domain into. */
        const int n_in;
        const int n;
        const int n_full;
        quadrat * q;
        /** The grid spacing. */
        double h;
        /** The Neumann condition to apply at x=2. */
        double g;

        double* const node_pos;
        double* const omega;
        double* const omega_1;
        double* x_full;
        double* b_full;
        double ** a_vals;
        double ** a_short;
        int ** a_ind;
        double ** bounds;

        cubic_1d_alt_C2(int n_);
        ~cubic_1d_alt_C2();

        void assemble_b();
        void assemble_a();

        void write_out(char prefix[], int N_test);

        double phi_C2(double x_in, int i);
        double grad_phi_C2(double x_in, int i);

    private:
        virtual void mul_A(double *in,double *out);
        virtual void M_inv(double *in,double *out);
        double f_source(double xx);
        double min(double R1, double R2)
        {
          if (R1 < R2) {return R1;}
          else {return R2;}
        }
        double max(double R1, double R2)
        {
          if (R1 < R2) {return R2;}
          else {return R1;}
        }

};

#endif
