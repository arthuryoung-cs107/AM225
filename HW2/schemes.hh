#ifndef SCHEMES_HH  /* Include guard */
#define SCHEMES_HH

#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "gsl/gsl_rng.h"

#include "Cash_Karp.hh"
#include "Geng.hh"
#include "Cash_Karp_GSL.hh"

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

class Kuramoto2D : public Cash_Karp_GSL
{
public:
  gsl_rng * T1;
  double J, K;
  Kuramoto2D(double lambda_in1, double lambda_in2, double J_in, double K_in, int tag) : Cash_Karp_GSL(1250, 3)
  {
    T1 = gsl_rng_alloc(gsl_rng_taus);
    J = J_in;
    K = K_in;
    writeout_width = dims*dof + 1;
    int i, j;
    double r_val, theta_val;

    a_tol = lambda_in1;
    r_tol = lambda_in2;
    memset(prefix, 0, 99);
    memset(prefix_dense, 0, 149);
    memset(specfile, 0, 199);
    memset(specfile_dense, 0, 249);

    snprintf(prefix, 100, "./dat_dir/prob5_2DKuramoto_sim%d", tag );
    snprintf(specfile, 200, "%s.aydat", prefix);
    snprintf(prefix_dense, 150, "%s_dense_output", prefix);
    snprintf(specfile_dense, 250, "%s.aydat", prefix_dense);

    for ( i = 0; i < dof; i++)
    {
      gsl_matrix_set(y_init, i, 0, 2*M_PI*gsl_rng_uniform(T1) );
      theta_val = 2*M_PI*gsl_rng_uniform(T1);
      r_val = gsl_rng_uniform(T1);
      gsl_matrix_set(y_init, i, 1, r_val*cos(theta_val) );
      gsl_matrix_set(y_init, i, 2, r_val*sin(theta_val) );
    }

  }
  ~Kuramoto2D()
  {
    gsl_rng_free(T1);
  }
  double omega_eval(gsl_matrix * y_in, int i)
  {
    double diff1, diff2;
    double acc = 0;
    for (int j = 0; j < i; j++)
    {
      diff1 = gsl_matrix_get(y_in, j, 1) - gsl_matrix_get(y_in, i, 1);
      diff2 = gsl_matrix_get(y_in, j, 2) - gsl_matrix_get(y_in, i, 2);
      acc += sin(gsl_matrix_get(y_in, j, 0) - gsl_matrix_get(y_in, i, 0) )/(sqrt( diff1*diff1 + diff2*diff2 ));
    }
    for (int j = i + 1; j < dof; j++)
    {
      diff1 = gsl_matrix_get(y_in, j, 1) - gsl_matrix_get(y_in, i, 1);
      diff2 = gsl_matrix_get(y_in, j, 2) - gsl_matrix_get(y_in, i, 2);
      acc += sin(gsl_matrix_get(y_in, j, 0) - gsl_matrix_get(y_in, i, 0) )/(sqrt( diff1*diff1 + diff2*diff2 ));
    }
    acc *= (K/( (double) dof ));
    return acc;
  }
  void x_eval(gsl_matrix * y_out, gsl_matrix * y_in, int i)
  {
    int j;
    double diff1, diff2;
    double acc1 = 0;
    double acc2 = 0;
    for ( j = 0; j < i; j++)
    {
      diff1 = gsl_matrix_get(y_in, j, 1) - gsl_matrix_get(y_in, i, 1);
      diff2 = gsl_matrix_get(y_in, j, 2) - gsl_matrix_get(y_in, i, 2);
      acc1 += (diff1/(sqrt(diff1*diff1 + diff2*diff2)))*(1.0 + J*cos(gsl_matrix_get(y_in, j, 0) - gsl_matrix_get(y_in, i, 0)) ) - (diff1/(diff1*diff1 + diff2*diff2));
      acc2 += (diff2/(sqrt(diff1*diff1 + diff2*diff2)))*(1.0 + J*cos(gsl_matrix_get(y_in, j, 0) - gsl_matrix_get(y_in, i, 0)) ) - (diff2/(diff1*diff1 + diff2*diff2));
    }
    for ( j = i + 1; j < dof; j++)
    {
      diff1 = gsl_matrix_get(y_in, j, 1) - gsl_matrix_get(y_in, i, 1);
      diff2 = gsl_matrix_get(y_in, j, 2) - gsl_matrix_get(y_in, i, 2);
      acc1 += (diff1/(sqrt(diff1*diff1 + diff2*diff2)))*(1.0 + J*cos(gsl_matrix_get(y_in, j, 0) - gsl_matrix_get(y_in, i, 0)) ) - (diff1/(diff1*diff1 + diff2*diff2));
      acc2 += (diff2/(sqrt(diff1*diff1 + diff2*diff2)))*(1.0 + J*cos(gsl_matrix_get(y_in, j, 0) - gsl_matrix_get(y_in, i, 0)) ) - (diff2/(diff1*diff1 + diff2*diff2));
    }
    acc1 *= 1/((double) dof );
    acc2 *= 1/((double) dof );
    gsl_matrix_set(y_out, i, 1, acc1);
    gsl_matrix_set(y_out, i, 2, acc2);
  }
  virtual void eval(double time, gsl_matrix * y_in, gsl_matrix * y_out)
  {
    # pragma omp parrallel for
      for (int i = 0; i < dof; i++)
      {
        gsl_matrix_set(y_out, i, 0, omega_eval(y_in, i));
        x_eval(y_out, y_in, i);
      }
    evals++;
  }
  virtual void write_out()
  {
    int i, j;
    double val;
    fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
    for ( i = 0; i < dof; i++)
    {
      for ( j = 0; j < dims; j++)
      {
        val = gsl_matrix_get(y0, i, j);
        fwrite(&(val), sizeof(double), 1, out_file_ptr);
      }
    }
  }
  void solve(double t_start, double t_end, int max_steps)
  {
    double deltaT = (t_end-t_start)/((double) max_steps);
    int out  = dense_solve(t_start, t_end, max_steps);
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
      writeout_width = dof + 1;
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
    virtual void write_out()
    {
      fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
      for (int i = 0; i < dof; i++) fwrite(&(y_it[i]), sizeof(double), 1, out_file_ptr) ;
    }
    void solve(double t_start, double t_end)
    {
      solve_fixed(t_start, t_end, max_steps_Brus);
    }
};
class Galaxy_Geng : public Geng
{
  public:
    int out_flag;
    double A, Omega, C, a, b, c;
    double * y_old, *y_work2;
    double t_old;
    Galaxy_Geng(double * params_in, int tag) : Geng(6),
    A(params_in[0]), Omega(params_in[1]), C(params_in[2]),
    a(params_in[3]), b(params_in[4]), c(params_in[5]), y_old(new double[6]), y_work2(new double[6])
    {
      memset(prefix, 0, 99);
      memset(specfile, 0, 199);
      snprintf(prefix, 100, "./dat_dir/prob6_Galaxy_data_sim%d", tag);
      snprintf(specfile, 200, "%s.aydat", prefix);
      out_file_ptr = fopen(specfile, "wb");
      // position, then momentum
      y_init[0] = 2.5;
      y_init[1] = 0.0;
      y_init[2] = 0.0;
      y_init[3] = 0.0;
      y_init[5] = 0.2;

      y_init[4] = solve_p2();
      writeout_width = dof + 2;
    }
    ~Galaxy_Geng()
    {
      delete [] y_old;
      delete [] y_work2;
      fclose(out_file_ptr);
    }
    virtual void eval(double time_in, double * y_in, double * y_out)
    {
      double q1 = y_in[0];
      double q2 = y_in[1];
      double q3 = y_in[2];
      double p1 = y_in[3];
      double p2 = y_in[4];
      double p3 = y_in[5];

      y_out[0] = p1 + Omega*q2;
      y_out[1] = p2 - Omega*q1;
      y_out[2] = p3;
      y_out[3] = -( -Omega*p2 + A*( (2.0*q1)/(a*a) )/(C + (q1*q1)/(a*a) + (q2*q2)/(b*b) + (q3*q3)/(c*c) ) );
      y_out[4] = -(  Omega*p1 + A*( (2.0*q2)/(b*b) )/(C + (q1*q1)/(a*a) + (q2*q2)/(b*b) + (q3*q3)/(c*c) ) );
      y_out[5] = -(             A*( (2.0*q3)/(c*c) )/(C + (q1*q1)/(a*a) + (q2*q2)/(b*b) + (q3*q3)/(c*c) ) );

      eval_count++;
    }
    virtual void init()
    {
      for (int i = 0; i < dof; i++) y_it[i] = y_init[i];
    }
    virtual void write_out()
    {
      if (out_flag == 1)
      {
        double Ham_eval = Hamiltonian(y_it);
        fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
        fwrite(&(Ham_eval), sizeof(double), 1, out_file_ptr);
        for (int i = 0; i < dof; i++) fwrite(&(y_it[i]), sizeof(double), 1, out_file_ptr) ;
      }
      else
      {
        if ( (y_it[1] > 0) && (y_it[4] > 0) )
        {
          if ( (y_old[1] < 0 ) )
          {
            double y2diff = y_it[1] - y_old[1];
            double twork = (y_it[1]*(t_old-t_it)/y2diff) + t_it;
            for (int i = 0; i < dof; i++)
            {
              if (i == 1)
              {
                y_work2[i] = 0;
              }
              else
              {
                y_work2[i] = (y_it[1]*(y_old[i] - y_it[i])/y2diff) + y_it[i] ;
              }
            }

            double Ham_eval = Hamiltonian(y_work2);
            fwrite(&(twork), sizeof(double), 1, out_file_ptr);
            fwrite(&(Ham_eval), sizeof(double), 1, out_file_ptr);
            for (int i = 0; i < dof; i++) fwrite(&(y_work2[i]), sizeof(double), 1, out_file_ptr) ;
          }
        }
        t_old = t_it;
        for (int i = 0; i < dof; i++) y_old[i] = y_it[i] ;
      }
    }
    void solve(double t_end, int out_flag_in)
    {
      out_flag = out_flag_in;
      int steps_ = ( (int) t_end ) * ( (int) 20 );
      solve_fixed(0.0, t_end, steps_);
    }
    double solve_p2()
    {
      double k1, k2, k3;
      double q1 = y_init[0];
      double q2 = y_init[1];
      double q3 = y_init[2];
      double p1 = y_init[3];
      double p3 = y_init[5];

      k1 = 0.5;
      k2 = -Omega*q1;
      k3 = 0.5*(p1*p1 + p3*p3) + Omega*q2*p1 + grav_potential(q1, q2, q3) - 2.0;
      return (-k2 + sqrt( k2*k2 - 4*k1*k3 ) )/(2*k1);
    }
    double grav_potential(double q1, double q2, double q3)
    {
      return A*log(C + (q1*q1/(a*a) ) + (q2*q2/(b*b)) + (q3*q3/(c*c)) );
    }
    double Hamiltonian(double * y_in)
    {
      double q1 = y_in[0];
      double q2 = y_in[1];
      double q3 = y_in[2];
      double p1 = y_in[3];
      double p2 = y_in[4];
      double p3 = y_in[5];

      return 0.5*(p1*p1 + p2*p2 + p3*p3) + Omega*(p1*q2-p2*q1) + grav_potential(q1, q2, q3);
    }


};


#endif
