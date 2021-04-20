#ifndef RITZ_GALERK_HH
#define RITZ_GALERK_HH

#include "conj_grad.hh"

#include <cstdio>
#include <cmath>

/** Class for solving an elliptic PDE problem over the domain [1,2] using
 * piecewise cubic basis functions. */
class Ritz_Galerk : public conj_grad {
    public:
        /** The number of intervals to divide the domain into. */
        const int n;
        /** The source function. */
        double* const f;
        /** The grid spacing. */
        const double h;
        /** The Neumann condition to apply at x=2. */
        double g;
        Ritz_Galerk(int n_) : conj_grad(3*n_),
           n(n_), f(new double[3*n+1]), h(1./3/n) {}
        virtual ~Ritz_Galerk() {delete [] f;}
        void init_const();
        void init_slope();
        void init_mms();
        double l2_norm_mms();
        void print_matrix();
        void print(FILE *fp);
        void print(const char* filename);
        inline void print() {print(stdout);}
        virtual void mul_A(double *in,double *out);
    private:
        inline double mms_dsq(double xx,double s) {
            double del=exp(1-xx)*sin(5*M_PI*xx)-s;
            return del*del;
        }
        void assemble_b();
};

#endif
