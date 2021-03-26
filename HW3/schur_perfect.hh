#include <functional>
#include "poisson_fft_AY.hh"
#include "conj_grad.hh"

class schur_perfect : conj_grad
{
    public:
        schur(int n_);
        ~schur();
        void solve(const std::function<double(double,double)>& f);
        virtual void mul_A(double* in, double* out);
        void print_solution();
        /** The number of interior gridpoints per subdomain in one dimension */
        const int n;
        /** The total number of interior gridpoints per subdomain */
        const int nn;
        /** The number of gridpoints on the glue */
        const int ng;
        /** The grid spacing */
        const double h;
        /** The inverse grid spacing squared */
        const double ih2;
    private:
        /** The solver for subdomain 1 */
        poisson_fft grid1;
        /** The solver for subdomain 2 */
        poisson_fft grid2;
        /** The source term on subdomain 1 */
        double* f1;
        /** The source term on subdomain 2 */
        double* f2;
        /** The source term on the glue */
        double* fg;
};
