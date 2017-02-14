#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <extern_global.h>
#include <file_manipulation.h>
#include <linalg.h>
#include <read_ini_file.h>
#include <runtime_messages.h>
#include <vectors.h>
#include <numerics.h>
#include <statistics.h>

struct BiasData {
  int Nk;
  double *k;
  double *b;
  double *bsigma;
};

void InitialiseCole05Data(char parfile[], struct BiasData * data)
{
  char datafile[STD_STR];
  if(!read_string_from_file(parfile,"biasdatafile",datafile)) read_err("biasdatafile");
  if(!check_status(status)) mpi_exit(status);

  data->Nk = count_lines(datafile);
  data->k = allocate_double_vector(data->Nk);
  data->b = allocate_double_vector(data->Nk);
  data->bsigma = allocate_double_vector(data->Nk);
  int i;
  FILE * f = fopen(datafile,"r");
  for(i=0; i<data->Nk; i++)
  {
    fscanf(f,"%lf %lf %lf",&data->k[i],&data->b[i],&data->bsigma[i]);
  }
  fclose(f);
}

void QuitCole05(struct BiasData * data)
{
  free(data->k);
  free(data->b);
  free(data->bsigma);
}

double chi2_Cole05(double *params, int Nparams, struct BiasData * data)
{
  // Make parameters readable
  double b0,A,Q;
  b0 = params[0];
  A = params[1];
  Q = params[2];

  double Mod[data->Nk];
  double deviation[data->Nk];

  int i;
  double chi2 = 0.0;

  for (i=0; i<data->Nk; i++)
  {
    Mod[i] = b0 * sqrt( (1.0+Q*(data->k[i]*data->k[i])) / (1.0+A*data->k[i]) );
    deviation[i] = data->b[i]-Mod[i];
  }

  for(i=0; i<data->Nk; i++)
  {
    // for(j=0;j<N;j++)
    // {
      // chi2 += deviation[i] * D->InvC[i][j] * deviation[j];
      chi2 += deviation[i] * deviation[i] / (data->bsigma[i] * data->bsigma[i]) ;
    // }
  }

  if(!check_for_nan_and_inf(chi2,"Cole05 chi2")) return -1;

  #ifdef DEBUG
  printf("chi2 = %e\n",chi2);
  #endif

  return chi2;
}
