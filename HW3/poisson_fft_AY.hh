#ifndef POISSON_FFT_AY
#define POISSON_FFT_AY

#include <cstdio>
#include <cmath>
#include <string.h>
#include <functional>

#include <fftw3.h>

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
}

class poisson_fft {
    public:
        /** The number of gridpoints in one dimension. */
        const int n;
        /** The total number of gridpoints. */
        const int nn;
        /** The grid spacing. */
        const double h;
        /** The discretized source term in the Poisson equation. */
        double* const f;
        /** The discretized solution to the Poisson equation. */
        double* const v;
        /** The frequency domain. */
        double* const w;
        poisson_fft(int n_, double h_);
        ~poisson_fft();
        void init(const std::function<double(double,double)>& f_in);
        void solve();
        void init_mms();
        double l2_error_mms();
        void print(bool solution);
        void output_solution( char prefix[], double * v_in);

    private:
        /** An array holding the eigenvalues of the one-dimensional Poisson
         * matrix T_N. */
        double* const lam;
        /** The FFTW plan for converting the source term into the frequency
         * domain. */
        fftw_plan plan_fwd;
        /** The FFTW plan for converting the frequency domain back to the
         * solution. */
        fftw_plan plan_bck;

        void write_out(char prefix[], int n, double * M);

};

#endif
