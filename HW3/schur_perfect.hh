#include <functional>
#include "poisson_fft_AY.hh"
#include "conj_grad.hh"
#include "square_specs.hh"

class schur_perfect : conj_grad
{
    public:
        int n_sqrs, n_glue, nn_glue, N;
        int * n_vec, *nn_vec, **mat_indices, ** vec_indices ; 
        double h, ih2;
        poisson_fft ** grids;
        double * f_full, * v_sol, ** fk_full, *fk;
        double *** A_glue;

        schur_perfect(square_specs * S);
        ~schur_perfect();
        void solve(const std::function<double(double,double)>& f);
        virtual void mul_A(double* in, double* out);
        void print_solution();

    private:

};
