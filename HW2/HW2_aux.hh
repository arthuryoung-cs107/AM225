#ifndef HW2_AUX_HH  /* Include guard */
#define HW2_AUX_HH

#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "Cash_Karp.hh"

class Brusselator : public Cash_Karp
{
public:
  Brusselator(double lambda_in1, double lambda_in2, int tag) : Cash_Karp(2)
  {
    a_tol = lambda_in1;
    r_tol = lambda_in2;
    memset(prefix, 0, 99);
    memset(prefix_dense, 0, 149);
    memset(specfile, 0, 199);
    memset(specfile_dense, 0, 249);

    snprintf(prefix, 100, "./dat_dir/prob2_Bruss_results_lambda%d", tag );
    snprintf(specfile, 200, "%s.aydat", prefix);
    snprintf(prefix_dense, 150, "%s_dense_output", prefix);
    snprintf(specfile_dense, 250, "%s.aydat", prefix_dense);
    out_file_ptr = fopen(specfile, "wb");
    out_file_dense_ptr = fopen(specfile_dense, "wb");

    y_init[0] = 1.5;
    y_init[1] = 3.0;
  }
  ~Brusselator()
  {
    fclose(out_file_dense_ptr);
    fclose(out_file_ptr);
  }
  virtual void eval(double time, double * y_in, double * y_out)
  {
    y_out[0] = 1 + y_in[0]*(y_in[0]*y_in[1]-4.0);
    y_out[1] = y_in[0]*(3.0 - y_in[0]*y_in[1] );
    evals++;
  }

};

class Two_Comp : public Cash_Karp
{
public:
  Two_Comp(double lambda_in1, double lambda_in2, int tag) : Cash_Karp(2)
  {
    a_tol = lambda_in1;
    r_tol = lambda_in2;
    memset(prefix, 0, 99);
    memset(prefix_dense, 0, 149);
    memset(specfile, 0, 199);
    memset(specfile_dense, 0, 249);

    snprintf(prefix, 100, "./dat_dir/prob2_TwoComp_results_lambda%d", tag );
    snprintf(specfile, 200, "%s.aydat", prefix);
    snprintf(prefix_dense, 150, "%s_dense_output", prefix);
    snprintf(specfile_dense, 250, "%s.aydat", prefix_dense);
    out_file_ptr = fopen(specfile, "wb");
    out_file_dense_ptr = fopen(specfile_dense, "wb");

    y_init[0] = 1.0;
    y_init[1] = 0.0;
  }
  ~Two_Comp()
  {
    fclose(out_file_dense_ptr);
    fclose(out_file_ptr);
  }
  virtual void eval(double time_in, double * y_in, double * y_out)
  {
    y_out[0] = -time_in*y_in[1];
    y_out[1] = time_in*y_in[0];
    evals++;
  }

};


#endif
