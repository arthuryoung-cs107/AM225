#ifndef CUBIC_1D_ALT_HH
#define CUBIC_1D_ALT_HH

#include "conj_grad.hh"
#include "quadrat.hh"

#include <cstdio>
#include <cmath>

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
        double a_C, a_J;

        double* const node_pos;
        double* const omega;

        cubic_1d_alt(int n_);
        ~cubic_1d_alt();

        void assemble_b();
        virtual void mul_A(double *in,double *out);
    private:
        double phi_C1(double x_in, int i);
        double f_source(double xx);
        double grad_phi_C1(double x_in, int i);

};

#endif
