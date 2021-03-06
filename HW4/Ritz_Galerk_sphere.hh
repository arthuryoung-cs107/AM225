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
#include "quadrat.hh"

#include <cstdio>
#include <cmath>

class Ritz_Galerk_sphere : public conj_grad {
    public:

        const int n, dof, n_full, n_q;

        double * xy, * vw, ** xy_coords, ** vw_coords;

        const double h;

        Ritz_Galerk_sphere(int n_);
        ~Ritz_Galerk_sphere();

        virtual void mul_A(double *in,double *out);
        void assemble_b(const std::function<double(double,double)>& f_source);
        void write_out(char prefix_x[], char prefix_y[], char prefix_v[], char prefix_w[], char prefix_u[], int N);
        void map(double * x_in, double * v_out);
        double phi_eval(double x_loc, double y_loc);

    private:
        AYmat * Jac, * a_vals, * Jac_inv, * dphi_i, * dphi_j, * work1, * work2;
        int * a_count, ** a_indices;
        quadrat * q;

        void Jac_eval(double * x_in, AYmat * mat_out );
        double Jac_det(AYmat * Jac_in);
        void grad_phi_eval(double x_loc, double y_loc, AYmat * grad_out);
        void assemble_a();
        double a_prod(int i, int j);
        void Jac_invert( AYmat * Jac_in, AYmat * Jac_inverted );
};

#endif
