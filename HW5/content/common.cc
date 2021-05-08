#include <cstdlib>

#include "common.hh"

#define NR_END 1
#define FREE_ARG char*

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* safe_fopen(const char* filename,const char* mode) {
    FILE *temp=fopen(filename,mode);
    if(temp==NULL) fprintf(stderr,"Error opening file \"%s\"",filename);
    return temp;
}

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void fatal_error(const char *p,int code) {
    fprintf(stderr,"Error: %s\n",p);
    exit(code);
}

void aysml_gen(char name[], int m, int n)
{
  char specfile[300];
  memset(specfile, 0, 299);
  snprintf(specfile, 300, "%s.aysml", name);
  FILE * aysml_file = fopen(specfile, "w");
  fprintf(aysml_file, "%d %d", m, n);
  fclose(aysml_file);
}

void fprintf_matrix(double ** matrix, int M, int N, char name[])
{
  int i, j;
  char specfile[200];
  memset(specfile, 0, 199);
  snprintf(specfile, 200, "%s.aydat", name);
  FILE * data_file = fopen(specfile, "wb");

  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      fwrite(&((*matrix)[i*N + j]), sizeof(double), 1, data_file);
    }
  }
  fclose(data_file);

  aysml_gen( name, M, N);

}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
