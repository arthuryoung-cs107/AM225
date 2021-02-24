#ifndef Cash_Karp_HH  /* Include guard */
#define Cash_Karp_HH

class Cash_Karp
{
public:
  const int dof;
  const double facmax, facmin, fac;
  double atol, rotol; 
  double *k1, *k2, *k3, *k4, *k5, *k6;
  double *w, *w_it;

  Cash_Karp(int dof_in);
  ~Cash_Karp();
  int solve(double t_start, double t_end, double ** outputs);
  inline int solve(double t_start, double t_end)
  {
    return solve(t_start, t_end, NULL);
  }
  virtual void eval(double time, double * y_in, double * y_out) = 0;
}

#endif
