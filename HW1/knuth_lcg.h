#ifndef knuth_lcg_H  /* Include guard */
#define knuth_lcg_H

u_int64_t lcg_uni(u_int64_t *lcg_carry); // call this evertime you want a random integer, or rand()
u_int64_t lcg_fwd(u_int64_t seed,u_int64_t jump); // eq to srand(seed)
double lcg_sze(); // replaces randmax
double boxmuller_knuth(double mean, double variance, u_int64_t * carry, double sze);
double uniform_random(double low, double high, u_int64_t * carry);

#endif // knuth_lcg_H
