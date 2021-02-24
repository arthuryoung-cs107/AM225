#ifndef HW2_AUX_HH  /* Include guard */
#define HW2_AUX_HH

class Brusselator : public Cash_Karp
{
public:
  brusselator() : Cash_Karp(2)
  {
    w[0] = 1.5;
    w[1] = 3.0;
  }
  virtual void eval(double time, double * y_in, double * y_out)
  {
    // time doesnt matter for brusselator, implicitly tied in
    double &y1=*y_in,&y2=y_in[1];
    y_out[0]=1+y1*(y1*y2-4);
    y_out[1]=y1*(3-y1*y2);
  }

}



#endif
