#include "Geng.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>

extern "C"
{
  #include "auxiliary_functions.h"
}


Geng::Geng(int dof_in) :  dof(dof_in),
y_init(new double[dof]),
y_it(new double[dof]),
y_final(new double[dof]),
y_work1(new double[dof]),
k1(new double[dof]), k2(new double[dof]), k3(new double[dof]),
k1b(new double[dof]), k2b(new double[dof]), k3b(new double[dof]) {}

Geng::~Geng()
{
  delete [] y_init;
  delete [] y_it;
  delete [] y_final;
  delete [] y_work1;

  delete [] k1;
  delete [] k2;
  delete [] k3;
  delete [] k1b;
  delete [] k2b;
  delete [] k3b;
}


void Geng::solve_fixed(double t_start, double t_end, int max_steps)
{
  double del_t = ( t_end - t_start )/( (double) max_steps );
  t_it = t_start;
  init();
  write_out();
  step_count = 1;

  for ( step_count = 2; step_count <= max_steps; step_count++)
  {
    step(del_t);
    write_out();
  }
  aysml_gen(prefix, max_steps, dof+1);
}

void Geng::step(double del_t)
{
  int i, it;
  const double root6 = sqrt(6.0);
  double delsqr, del;
  double * switch_ptr;

  const double c1 = ( 4.0 - root6 )/( 10.0 ),
               c2 = ( 4.0 + root6 )/( 10.0 ),
               c3 = 1.0,

               b1 = ( 16.0 - root6)/( 36.0 ),
               b2 = ( 16.0 + root6)/( 36.0 ),
               b3 = ( 1.0 )/( 9.0 ),

               a11 = (16.0 - root6 )/( 72.0 ),
               a12 = ( 328.0 - 167.0*root6 )/( 1800.0 ),
               a13 = ( -2.0 + 3.0*root6 )/( 450.0 ),
               a21 = ( 328.0 + 167.0*root6 )/( 1800.0 ),
               a22 = ( 16.0 + root6 )/( 72.0 ),
               a23 = ( -2.0 - 3.0*root6 )/( 450.0 ),
               a31 = ( 85.0 - 10.0*root6 )/( 180.0 ),
               a32 = ( 85.0 + 10.0*root6 )/( 180.0 ),
               a33 = ( 1.0 )/( 18.0 );

   for ( i = 0; i < dof; i++) k1[i]=k2[i]=k3[i] = 0.0 ;
   it = 0;
   do // fixed point iteration
   {
     if(++it>1000)
     {
       printf("implicit Runge-Kutta scheme did not converge. Halp\n");
       exit(1);
     }

     for ( i = 0; i < dof; i++) y_work1[i] = y_it[i] + del_t*(a11*k1[i] + a12*k2[i] + a13*k3[i]) ;
     eval(t_it + del_t*c1, y_work1, k1b);

     for ( i = 0; i < dof; i++) y_work1[i] = y_it[i] + del_t*(a21*k1[i] + a22*k2[i] + a23*k3[i]) ;
     eval(t_it + del_t*c2, y_work1, k2b);

     for ( i = 0; i < dof; i++) y_work1[i] = y_it[i] + del_t*(a31*k1[i] + a32*k2[i] + a33*k3[i]) ;
     eval(t_it + del_t*c3, y_work1, k3b);


     delsqr = 0;
     for ( i = 0; i < dof; i++) // l2 norm
     {
       del = k1[i] - k1b[i];
       delsqr += del*del;

       del = k2[i] - k2b[i];
       delsqr += del*del;

       del = k3[i] - k3b[i];
       delsqr += del*del;
     }

     switch_ptr = k1b;
     k1b = k1;
     k1 = switch_ptr;

     switch_ptr = k2b;
     k2b = k2;
     k2 = switch_ptr;

     switch_ptr = k3b;
     k3b = k3;
     k3 = switch_ptr;

   } while(delsqr > 1e-15);

   // Following convergence, update y

   for ( i = 0; i < dof; i++) y_it[i] += del_t*( b1*k1[i] + b2*k2[i] + b3*k3[i] );
   t_it += del_t;

}

void Geng::write_out()
{
  fwrite(&(t_it), sizeof(double), 1, out_file_ptr);
  for (int i = 0; i < dof; i++) fwrite(&(y_it[i]), sizeof(double), 1, out_file_ptr) ;
}
