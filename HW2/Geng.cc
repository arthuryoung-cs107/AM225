#include "Geng.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>

Geng::Geng(int dof_in) :  dof(dof_in),
y_final(new double[dof]),
k1(new double[dof]), k2(new double[dof]), k3(new double[dof]),
k1b(new double[dof]), k2b(new double[dof]), k3b(new double[dof]),
k1c(new double[dof]), k2b(new double[dof]), k3c(new double[dof]) {}

Geng::~Geng()
{
  delete [] y_final;
  delete [] k1;
  delete [] k2;
  delete [] k3;
  delete [] k1b;
  delete [] k2b;
  delete [] k3b;
  delete [] k1c;
  delete [] k2c;
  delete [] k3c;
}

void Geng::solve_fixed(double t_start, double t_end, double del_t)
{
  init();

  for (int i = 0; i < count; i++)
  {
    step(dt);
  }
}

void Geng::step(double del_t)
{
  
}
