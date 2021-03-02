#ifndef HW2_AUX_HH  /* Include guard */
#define HW2_AUX_HH

#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "Cash_Karp.hh"

class Brusselator : public Cash_Karp
{
public:
  Brusselator() : Cash_Karp(2)
  {
    memset(prefix, 0, 99);
    memset(specfile, 0, 199);
    snprintf(prefix, 100, "./dat_dir/prob2_Bruss_results");
    snprintf(specfile, 200, "%s.aydat", prefix);
    out_file_ptr = fopen(specfile, "wb");

    y_init[0] = 1.5;
    y_init[1] = 3.0;
  }
  virtual void eval(double time, double * y_in, double * y_out)
  {
    y_out[0] = 1 + y_in[0]*(y_in[0]*y_in[1]-4.0);
    y_out[1] = y_in[0]*(3.0 - y_in[0]*y_in[1] );
  }

};



#endif
