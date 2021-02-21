#ifndef auxiliary_functions_H  /* Include guard */
#define auxiliary_functions_H

void aysml_gen(char name[], int m, int n);
void fprintf_matrix(double ** matrix, int M, int N, char prefix[]);

// access functions
float t_acc(float * ptr, int i);
float x_acc(float * ptr, int i);
float xp_acc(float * ptr, int i);
double t_acc_d(double * ptr, int i);
double x_acc_d(double * ptr, int i);
double xp_acc_d(double * ptr, int i);
void t_stor(float * ptr, int i, float t);
void x_stor(float * ptr, int i, float x);
void xp_stor(float * ptr, int i, float xp);
void t_stor_d(double * ptr, int i, double t);
void x_stor_d(double * ptr, int i, double x);
void xp_stor_d(double * ptr, int i, double xp);

// math functions
double dmax_element(double * vector, int vector_low, int vector_high );
double dmin_element(double * vector, int vector_low, int vector_high ); //assumes index zero
int imax_element(int * ivector, int vector_low, int vector_high );
int imin_element(int * ivector, int vector_low, int vector_high );
void zerom_init(double **a, int arl, int arh, int acl, int ach);
void zeromint_init(int **a, int arl, int arh, int acl, int ach);
void onev_init(double *a, int arl, int arh);
void zerov_init(double *a, int arl, int arh);
void zerovint_init(int *a, int arl, int arh);
void dmatrix_mult(double ** a, int arl, int arh, int acl, int ach, double ** b, int brl, int brh, int bcl, int bch, double ** c );
void dmv_mult(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh, double * c);
void dmatrix_add(double ** a, int arl, int arh, int acl, int ach, double **b, int brl, int brh, int bcl, int bch, int sign, double ** result);
void dvector_add(double * a, int arl, int arh, double * b, int sign, double * result);
void dv_scalmult(double * a, int arl, int arh, double scalar, double * result);
void dm_scalmult(double ** a, int arl, int arh, int acl, int ach, double scalar, double ** result);
void mat2vec(double ** G, int nrl, int nrh, int ncl, int nch, double * g, int mrl, int mrh);
void vec2mat(double *g, int nrl, int nrh, double ** G, int mrl, int mrh, int mcl, int mch);
double norm_l2(double * input_vector, int nrl, int nrh);
double norm_frob(double ** X, int nrl, int nrh, int ncl, int nch);
double up_order(double a);


#endif // auxiliary_functions_H
