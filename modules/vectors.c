#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include <runtime_messages.h>
#include <extern_global.h>

int * allocate_int_vector(int nelems)
{
  int * vec = (int*)malloc(nelems*sizeof(int));
  if (vec==NULL) {frame("Bad memory allocation"); status = FAILURE;}
  return vec;
}

double * allocate_double_vector(int nelems)
{
  double * vec = (double*)malloc(nelems*sizeof(double));
  if (vec==NULL) {frame("Bad memory allocation"); status = FAILURE;}
  return vec;
}

char * allocate_char_vector(int nelems)
{
  char * vec = (char*)malloc(nelems*sizeof(char));
  if (vec==NULL) {frame("Bad memory allocation"); status = FAILURE;}
  return vec;
}

int ** allocate_int_matrix(int nrows, int ncols)
{
  int **m = (int **)malloc(nrows*sizeof(int*));
  int i;
  if (m!=NULL)
  {
    for (i=0; i<nrows; i++) m[i] = allocate_int_vector(ncols);
    if(!status) return NULL;
  }
  else
  {
    frame("Bad memory allocation"); status = FAILURE;
  }
  return m;
}

double ** allocate_double_matrix(int nrows, int ncols)
{
  double **m = (double **)malloc(nrows*sizeof(double*));
  int i;
  if (m!=NULL)
  {
    for (i=0; i<nrows; i++) m[i] = allocate_double_vector(ncols);
    if (!status) return NULL;
  }
  else
  {
    frame("Bad memory allocation"); status = FAILURE;
  }
  return m;
}

char ** allocate_char_matrix(int nrows, int ncols)
{
  char **m = (char **)malloc(nrows*sizeof(char*));
  int i;
  if (m!=NULL)
  {
    for (i=0; i<nrows; i++) m[i] = allocate_char_vector(ncols);
    if (!status) return NULL;
  }
  else
  {
    frame("Bad memory allocation"); status = FAILURE;
  }
  return m;
}

void deallocate_int_vector(int *vec)
{
  if(vec!=NULL) free(vec);
}

void deallocate_double_vector(double *vec)
{
  if(vec!=NULL) free(vec);
}

void deallocate_char_vector(char *vec)
{
  if(vec!=NULL) free(vec);
}

void deallocate_int_matrix(int **m, int nrows, int ncols)
{
  if(m!=NULL)
  {
    int i; for(i=0;i<nrows;i++) deallocate_int_vector(m[i]);
    free(m);
  }
}

void deallocate_double_matrix(double **m, int nrows, int ncols)
{
  if(m!=NULL)
  {
    int i; for(i=0;i<nrows;i++) deallocate_double_vector(m[i]);
    free(m);
  }
}

void deallocate_char_matrix(char **m, int nrows, int ncols)
{
  if(m!=NULL)
  {
    int i; for(i=0;i<nrows;i++) deallocate_char_vector(m[i]);
    free(m);
  }
}
