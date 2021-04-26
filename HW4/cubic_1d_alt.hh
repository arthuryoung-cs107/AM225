#ifndef CUBIC_1D_ALT_HH
#define CUBIC_1D_ALT_HH

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
class cubic_1d_alt : public conj_grad {
    public:
        /** The number of intervals to divide the domain into. */
        const int n;
        quadrat * q;
        /** The grid spacing. */
        const double h;
        /** The Neumann condition to apply at x=2. */
        double g;

        double* const node_pos;
        double* const omega;
        double ** a_vals;
        int * a_count;

        cubic_1d_alt(int n_);
        ~cubic_1d_alt();

        void assemble_b();
        void assemble_a();
        void write_out(char prefix[], int N_test);
        double phi_C1(double x_in, int i);
        double grad_phi_C1(double x_in, int i);

    private:
        virtual void mul_A(double *in,double *out);
        double f_source(double xx);

};

#endif
