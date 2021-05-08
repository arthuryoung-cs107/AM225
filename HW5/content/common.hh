#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>
#include <cmath>
#include <cstring>


FILE* safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);
void aysml_gen(char name[], int m, int n);
void fprintf_matrix(double ** matrix, int M, int N, char prefix[]);

double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

#endif
