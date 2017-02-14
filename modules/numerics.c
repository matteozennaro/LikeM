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
#include <runtime_messages.h>

double min_double_vec(int N, double *vec)
{
  double min = vec[0];
  int i;
  for (i=0;i<N;i++)
  {
    if(vec[i]<min) min=vec[i];
  }
  return min;
}

double max_double_vec(int N, double *vec)
{
  double max = vec[0];
  int i;
  for (i=0;i<N;i++)
  {
    if(vec[i]>max) max=vec[i];
  }
  return max;
}


double lin_interp(double x0, double *x, double *y, int dimension)
{
  if (x0 < x[0])
  {
    int i = 1;
    return ((y[i]-y[i-1])/(x[i]-x[i-1])*x0 + y[i] - (y[i]-y[i-1])/(x[i]-x[i-1])*x[i]);
  }
  else if (x0 > x[dimension-1])
  {
    int i = dimension-1;
    return ((y[i]-y[i-1])/(x[i]-x[i-1])*x0 + y[i] - (y[i]-y[i-1])/(x[i]-x[i-1])*x[i]);
  }
  else
  {
    int i=0;
    while (x0 > x[i]) i++;
    return ((y[i]-y[i-1])/(x[i]-x[i-1])*x0 + y[i] - (y[i]-y[i-1])/(x[i]-x[i-1])*x[i]);
  }
}

int check_for_nan_and_inf(double val, char whatsthis[])
{
  if(isnan(val))
	{
		char error[LONG_STR];
		sprintf(error,"Your %s is NaN!",whatsthis);
		frame(error);
		status = NUM_FAILURE;
		return 0;
	}

	if(isinf(val))
	{
		char error[LONG_STR];
		sprintf(error,"Your %s is infinity!",whatsthis);
		frame(error);
		status = NUM_FAILURE;
		return 0;
	}

  return 1;
}
