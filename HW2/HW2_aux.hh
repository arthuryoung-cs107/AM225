#ifndef HW2_AUX_HH  /* Include guard */
#define HW2_AUX_HH

#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "Cash_Karp.hh"
#include "Geng.hh"

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

class Brusselator_Geng : public Geng
{
  public:
    int max_steps_Brus;
    Brusselator_Geng(int max_steps, int tag1, int tag2) : Geng(2), max_steps_Brus(max_steps)
    {
      memset(prefix, 0, 99);
      memset(specfile, 0, 199);
      snprintf(prefix, 100, "./dat_dir/prob6_BrussGeng_results_steps_%de%d", tag1, tag2 );
      snprintf(specfile, 200, "%s.aydat", prefix);
      out_file_ptr = fopen(specfile, "wb");
      y_init[0] = 1.5;
      y_init[1] = 3.0;
    }
    ~Brusselator_Geng()
    {
      fclose(out_file_ptr);
    }
    virtual void eval(double time_in, double * y_in, double * y_out)
    {
      y_out[0] = 1 + y_in[0]*(y_in[0]*y_in[1]-4.0);
      y_out[1] = y_in[0]*(3.0 - y_in[0]*y_in[1] );
      eval_count++;
    }
    virtual void init()
    {
      for (int i = 0; i < dof; i++) y_it[i] = y_init[i];
    }
    void solve(double t_start, double t_end)
    {
      solve_fixed(t_start, t_end, max_steps_Brus);
    }
};
#endif
