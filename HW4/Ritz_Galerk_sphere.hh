#ifndef RITZ_GALERK_SPHERE_HH
#define RITZ_GALERK_SPHERE_HH

#include <cstdio>
#include <cstdlib>
// #include <cmath>
// #include <stdint.h>
// #include <string.h>
#include <functional>

#include "conj_grad.hh"
#include "AYmat.hh"

#include <cstdio>
#include <cmath>

/** Class for solving an elliptic PDE problem over the domain [1,2] using
 * piecewise cubic basis functions. */
class Ritz_Galerk_sphere : public conj_grad {
    public:
        /** The number of intervals to divide the domain into. */
        const int n, dof, n_full;
        /** The source function. */
        double * f;
        double * xy, * vw;
        double ** xy_coords;
        double ** vw_coords;
        /** The grid spacing. */
        const double h;

        Ritz_Galerk_sphere(int n_);
        ~Ritz_Galerk_sphere();

        virtual void mul_A(double *in,double *out);
    private:
        AYmat * Jac;

        void map(double * x_in, double * v_out);
        void Jac_eval(double * x_in, AYmat * mat_out );
        double Jac_det(AYmat * Jac_in);
        double phi_eval(double x_loc, double y_loc);
        void assemble_b(const std::function<double(double,double)>& f_source);
};

#endif
