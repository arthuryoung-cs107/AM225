#include "poisson_fft_AY.hh"

/** Initializes the class for solving the 2D Poisson problem on a square
 * [0,1]^2 using the fast Fourier transform.
 * \param[in] n the number of non-zero gridpoints in one direction. */
poisson_fft::poisson_fft(int n_, double h_)
    : n(n_), nn(n*n), h(h_), f(fftw_alloc_real(nn)),
    v(fftw_alloc_real(nn)), w(fftw_alloc_real(nn)), lam(new double[n]),
    plan_fwd(fftw_plan_r2r_2d(n,n,f,w,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE)),
    plan_bck(fftw_plan_r2r_2d(n,n,w,v,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE)) {

    // Initialize the table of eigenvalues
    const double fac=M_PI/(n+1);
    for(int j=0;j<n;j++) lam[j]=2*(1-cos(fac*(j+1)));
}

/** The class destructor frees the dynamically allocated memory, including the
 * FFTW plans, FFTW arrays, and eigenvalue table. */
poisson_fft::~poisson_fft() {
    fftw_destroy_plan(plan_bck);
    fftw_destroy_plan(plan_fwd);
    delete [] lam;
    fftw_free(w);
    fftw_free(v);
    fftw_free(f);
}

void poisson_fft::init(const std::function<double(double,double)>& f_in) {
    for(int j=0;j<n;j++)
    {
        double y_loc=( (double) j+1)*h;
        for(int i=0;i<n;i++) {
            double x_loc=( (double) i+1)*h;
            f[i+n*j] = f_in(x_loc, y_loc);
        }
    }
}

/** Solves the linear system using the fast Fourier transform. */
void poisson_fft::solve() {
    const double nor=1./(2*(n+1)),fac=h*h*nor*nor;

    // Perform the discrete sine transform using FFTW to convert the source
    // term into the frequency domain
    fftw_execute(plan_fwd);

    // Multiply each mode component by the corresponding scaling factor. Note
    // that a factor of nor^2 is included to deal with the overall scaling of
    // the FFTW routines.
    for(int j=0;j<n;j++) for(int i=0;i<n;i++) {
        w[i+n*j]*=fac/(lam[i]+lam[j]);
    }

    // Perform the discrete sine transform again to obtain the solution
    fftw_execute(plan_bck);
}

/** Prints either the source term or the solution as a grid of text values.
 * \param[in] solution whether to print the solution */
void poisson_fft::print(bool solution) {
    double *ptr=solution?v:f;
    for(int j=0;j<n;j++) {
        printf("%g",v[j]);
        for(int i=1;i<n;i++) printf(" %g",ptr[j+i*n]);
        putchar('\n');
    }
}

/** Outputs the solution in the Gnuplot 2D matrix format, padding the grid with
 * zeros to represent the Dirichlet boundary condition.
 * \param[in] filename the name of the file to write to. */
void poisson_fft::output_solution( char prefix[], double * v_in) {
    const int ne=n+2;
    double *fld=new double[(n+2)*(n+2)],*f2=fld;

    // Set first row to zero
    while(f2<fld+ne) *(f2++)=0.;

    // Copy field contents into output array, padding the start and end entries
    // with zeros
    for(int j=0;j<n;j++) {
        *(f2++)=0;
        for(int i=0;i<n;i++) *(f2++)=v_in[i+j*n];
        *(f2++)=0;
    }

    // Set last row to zero, call output routine, and free output array
    while(f2<fld+ne*ne) *(f2++)=0.;    write_out(prefix, ne, fld);
    delete [] fld;
}
void poisson_fft::write_out( char prefix[], int n, double * M)
{
  int i, j;
  FILE * out_file_ptr;
  char specfile[300];

  memset(specfile, 0, 299);
  snprintf(specfile, 200, "%s.aydat", prefix);
  out_file_ptr = fopen(specfile, "wb");

  for ( j = 0; j < n; j++)
  {
    for ( i = 0; i < n; i++) // walk
    {
      fwrite(&(M[i+n*j]), sizeof(double), 1, out_file_ptr);
    }
  }
  aysml_gen(prefix, n, n);

  fclose(out_file_ptr);
}
