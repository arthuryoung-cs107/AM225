#ifndef AYMAT_OPS_HH
#define AYMAT_OPS_HH

#include <cstdio>

#include "AYmat.hh"

void AYmat_mul(AYmat * A1, AYmat * A2, AYmat * A3);
void AYmat_mul_Strass(AYmat * A1, AYmat * A2, AYmat * A3);
void Strass_recurse(double ** A1, double ** A2, double ** A3, int N);

#endif
