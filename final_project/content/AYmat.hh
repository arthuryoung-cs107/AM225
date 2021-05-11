#ifndef AYMAT_HH
#define AYMAT_HH

#include <cstring>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

class AYmat {
    public:
      int M, N;
      int GSL_flag;
      double * A_ptr;
      gsl_matrix * A_gsl;

      AYmat(int M_, int N_);
      ~AYmat();

      void GSL_init();
      void GSL_send();

      void print_mat();
      void init_123();
      void init_0();
      void init_randuni();
      void set(int i, int j, double val);
      double get(int i, int j);

      AYmat * copy_gen();
      void copy_set(AYmat * X);

      AYmat * transpose_gen();
      void transpose_set(AYmat * XT);

      void add(AYmat * B, AYmat * C, double alpha, double beta);

      AYmat * mult_gen(AYmat * B, double alpha);
      void mult_set(AYmat * B, AYmat * C, double alpha, double beta);

      double inner(AYmat * B);
      double norm_frob();
      double norm_1();

      void svd(gsl_vector * S, gsl_matrix * V, gsl_vector * work);

    private:
      double ** AT;
};

AYmat * aysml_read(char name[]);


#endif
