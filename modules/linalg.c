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

#include <extern_global.h>
#include <read_ini_file.h>

void multiply_squared_matrices(double **M1, double **M2, double **RES, int N)
{
  int i,j,k;
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      RES[i][j] = 0.0;
      for(k=0;k<N;k++)
      {
        RES[i][j] += M1[i][k]*M2[k][j];
      }
    }
  }
}

void LineariseGSLMatrix(gsl_matrix * m, double * vec, int Nrows, int Ncols)
{
  int i,j;
  for(i=0;i<Nrows;i++)
  {
    for(j=0;j<Ncols;j++)
    {
      vec[i*Ncols+j] = gsl_matrix_get(m,i,j);
    }
  }
}
