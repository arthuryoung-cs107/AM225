#include <functional>
#include "poisson_fft_AY.hh"
#include "conj_grad.hh"
#include "square_specs.hh"

class schur_perfect : conj_grad
{
    public:
        int n_sqrs, n_glue, nn_glue, N;
        int * n_vec, *nn_vec, **mat_indices, ** vec_indices, ** coords ;
        double h, ih2;
        poisson_fft ** grids;
        double * f_full, * v_sol, * v_sol2,  ** fk_full, *fk, *Minv_mat, *A_check;
        double *** A_glue, ***A_glueT;

        schur_perfect(square_specs * S);
        ~schur_perfect();
        void init(const std::function<double(double,double)>& f);
        void solve_S();
        virtual void mul_A(double* in, double* out);
        virtual void M_inv(double *in,double *out);
        void output_solution( char prefix[], double * v_in);
        void write_out( char prefix[], int n, double * M);
    private:
        void init_preconditioning();
};
