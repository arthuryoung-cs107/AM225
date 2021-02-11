#ifndef HW1_AUX_HH
#define HW1_HEADER_HH

#include <cstdio>
#include <stdint.h>
#include "gsl/gsl_rng.h"

double random_uni(double low, double high, uint64_t * carry);
uint64_t casino_game(gsl_rng * T);
void cell_automaton(int m, int n, int number_threads);
double uniform_random(double low, double high, uint64_t * carry);

#endif
