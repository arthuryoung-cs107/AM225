#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

u_int64_t lcg_uni(u_int64_t *lcg_carry) // call this evertime you want a random integer, or rand()
{
  // **********************************************
  // * KNUTH's 64-BIT LCG RANDOM NUMBER GENERATOR *
  // **********************************************

  #ifdef KNUTH
    //ax+c mod m: m=2^64
    *lcg_carry= (__uint128_t) *lcg_carry*6364136223846793005ULL+1442695040888963407ULL; //knuth
  #else
    *lcg_carry=rand();
  #endif

  return *lcg_carry;
}

u_int64_t lcg_fwd(u_int64_t seed,u_int64_t jump) // eq to srand(seed)
{
  // *********************
  // * LCG JUMP FUNCTION *
  // *********************

  //being lazy here, should use the function  a^k x + (a^k-1)*c/(a-1)

  #ifdef KNUTH
    u_int64_t ix;
    u_int64_t ret;
    ret=0;
    for(ix=0;ix<jump;ix++)
    {
      ret=lcg_uni(&seed);
    }
    return ret;
  #else
    srand(seed);
    return seed;
  #endif
}

double lcg_sze() // replaces randmax
{
  // *******************************
  // * THE DOMAIN OF RANDOM NUMERS *
  // *******************************/

  #ifdef KNUTH
    return pow(2,64);
  #else
    return (double) RAND_MAX + 1.0;
  #endif
}
double boxmuller_knuth(double mean, double variance, u_int64_t * carry, double sze)
{
  double rand_uni1 = ((double) lcg_uni(carry))/sze; //numerator is double, forces double arithmatic
  double rand_uni2 = ((double) lcg_uni(carry))/sze;
  double box1 = -2 * log(rand_uni1);
  double gauss = pow(box1, 0.5)* cos(2*M_PI*rand_uni2);
  double absg = fabs(gauss);

  if (absg > 4)
  {
    gauss = 4 * gauss/absg;
  }

  gauss = (gauss * pow(variance, 0.5)) + mean;
  return gauss;
}

double uniform_random(double low, double high, u_int64_t * carry)
{
  double rand_uni = ((double) lcg_uni(carry))/lcg_sze(); //numerator is double, forces double arithmatic
  rand_uni = rand_uni*(high - low) + ((high - low)/2 - 0.5);
  return rand_uni;
}
