#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

// #include "blas.h"
// #include "omp.h"
#include "fluid_2d.hh"

// extern "C"
// {
//   #include "nrutil.h"
//   #include "auxiliary_functions.h"
// }
const char fn[]="ftest.out";

void color_init1(fluid_2d * m2in)
{
  double work;
  int check;

  int m = m2in->m; int n = m2in->n; int ml = m2in->ml;
  double ay = m2in->ay; double ax = m2in->ax; double dy = m2in->dy; double dx = m2in->dx;
  field * fm =  m2in->fm;


// #pragma omp parallel for
      for(int j=0;j<n;j++) // generic nodes
      {
          double y=ay+dy*(j+0.5);
          field *fp=fm+ml*j;
          for(int i=0;i<m;i++)
          {
              double x=ax+dx*(i+0.5);
              work = ((x*6.0) + 6.0 );
              check = (int) work;
              work = ((y*6.0) + 6.0 );
              check += (int) work;
              if (check%2 == 0)
              {
                fp->R = 1.0;
                fp->G = 1.0;
                fp->B = 1.0;
              }
              else
              {
                fp->R = 0.2;
                fp->G = 0.4;
                fp->B = 0.9;
              }
              fp++;
          }
      }
      for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) // left and right ghost nodes
      {
          fp[-2].R = fp[-1].R = fp[0].R;
          fp[-2].G = fp[-1].G = fp[0].G;
          fp[-2].B = fp[-1].B = fp[0].B;

          fp[m+2].R = fp[m+1].R = fp[m].R;
          fp[m+2].G = fp[m+1].G = fp[m].G;
          fp[m+2].B = fp[m+1].B = fp[m].B;
      }

      const int tl=2*ml,g=n*ml;
      for(field *fp=fm-2,*fe=fm+m+2;fp<fe;fp++) // top and bottom ghost nodes
      {
        fp[-tl].R = fp[-ml].R = fp[0].R;
        fp[-tl].G = fp[-ml].G = fp[0].G;
        fp[-tl].B = fp[-ml].B = fp[0].B;

        fp[g].R = fp[g+ml].R = fp[g-ml].R;
        fp[g].G = fp[g+ml].G = fp[g-ml].G;
        fp[g].B = fp[g+ml].B = fp[g-ml].B;
     }
}

void color_init2(fluid_2d * m2in)
{
  double work;
  int check;

  int m = m2in->m; int n = m2in->n; int ml = m2in->ml;
  double ay = m2in->ay; double ax = m2in->ax; double dy = m2in->dy; double dx = m2in->dx;
  field * fm =  m2in->fm;


// #pragma omp parallel for
      for(int j=0;j<n;j++) // generic nodes
      {
          double y=ay+dy*(j+0.5);
          field *fp=fm+ml*j;
          for(int i=0;i<m;i++)
          {
              double x=ax+dx*(i+0.5);
              work = ((x*6.0) + 6.0 );
              check = (int) work;
              work = ((y*6.0) + 6.0 );
              check += (int) work;
              if (check%2 == 0)
              {

              }
              else
              {
                fp->R = 0.2;
                fp->G = 0.4;
                fp->B = 0.9;
              }
              fp++;
          }
      }
      for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) // left and right ghost nodes
      {
          fp[-2].R = fp[-1].R = fp[0].R;
          fp[-2].G = fp[-1].G = fp[0].G;
          fp[-2].B = fp[-1].B = fp[0].B;

          fp[m+2].R = fp[m+1].R = fp[m].R;
          fp[m+2].G = fp[m+1].G = fp[m].G;
          fp[m+2].B = fp[m+1].B = fp[m].B;
      }

      const int tl=2*ml,g=n*ml;
      for(field *fp=fm-2,*fe=fm+m+2;fp<fe;fp++) // top and bottom ghost nodes
      {
        fp[-tl].R = fp[-ml].R = fp[0].R;
        fp[-tl].G = fp[-ml].G = fp[0].G;
        fp[-tl].B = fp[-ml].B = fp[0].B;

        fp[g].R = fp[g+ml].R = fp[g-ml].R;
        fp[g].G = fp[g+ml].G = fp[g-ml].G;
        fp[g].B = fp[g+ml].B = fp[g-ml].B;
     }
}

void prob1_part_a()
{

  // char prefix[200];
  // memset(prefix, 0, 199);
  // snprintf(prefix, 100, "./dat_dir/prob1");

  // Create the output directory for storing the simulation frames
  mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

  // Specify which fields should be outputted. 1: horizontal velocity, 2:
  // vertical velocity, 4: pressure.
  unsigned int fflags=1|2|4;

  // Construct the simulation class, setting the number of gridpoints, the
  // periodicity, and physical constants
  // fluid_2d f2d(256,256,true,true,-1,1,-1,1,0.002,1.,fflags,fn);
  fluid_2d m2d(256,256,true,true,-1,1,-1,1,0.002,1.,fflags,fn);

  // Initialize the tracers, and set the timestep based on multiplying the
  // maximum allowable by a padding factor
  // f2d.initialize(512,0.6);
  m2d.initialize(512,0.6);
  color_init1(&m2d);

  // Run the simulation for a specified duration, outputting snapshots at
  // regular intervals
  m2d.solve(0.5,200);


}

int main()
{

  prob1_part_a();

  return 0;
}
