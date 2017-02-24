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

struct DataPk {
  int Nk;
  double *k;
  double *pk;
  double *err;
};

void InitialiseCole05DataPk(char parfile[], struct DataPk * data, struct DataPk * lin)
{
  char datafile[STD_STR];
  if(!read_string_from_file(parfile,"pkdatafile",datafile)) read_err("pkdatafile");
  if(!check_status(status)) mpi_exit(status);

  data->Nk = count_lines(datafile)-1;
  data->k = allocate_double_vector(data->Nk);
  data->pk = allocate_double_vector(data->Nk);
  data->err = allocate_double_vector(data->Nk);
  int i;
  FILE * f = fopen(datafile,"r");
  char buf[2000];
  double val, psn;
  fscanf(f,"%lf %lf",&val,&psn);
  for(i=0; i<data->Nk; i++)
  {
    fscanf(f,"%lf %lf %lf %lf",&val,&data->k[i],&data->pk[i],&data->err[i]);
    data->err[i] = sqrt(2.0/data->err[i])*data->pk[i];
    data->pk[i] -= psn;
  }
  fclose(f);

  if(!read_string_from_file(parfile,"pklinfile",datafile)) read_err("pklinfile");
  if(!check_status(status)) mpi_exit(status);
  int nl,nh;
  nh = count_header_lines(datafile);
  nl = count_lines(datafile)-nh;
  double ktmp[nl],pktmp[nl];
  f = fopen(datafile,"r");
  for(i=0; i<nh; i++) fgets(buf,sizeof(buf),f);
  for(i=0; i<nl; i++)
  {
    fscanf(f,"%lf %lf",&ktmp[i],&pktmp[i]);
  }
  fclose(f);
  lin->Nk = data->Nk;
  lin->k = allocate_double_vector(lin->Nk);
  lin->pk = allocate_double_vector(lin->Nk);
  lin->err = allocate_double_vector(lin->Nk);
  for(i=0;i<lin->Nk;i++)
  {
    lin->k[i] = data->k[i];
    lin->pk[i] = lin_interp(lin->k[i],ktmp,pktmp,nl) / ((2.0*M_PI)*(2.0*M_PI)*(2.0*M_PI));
    lin->err[i] = 0.0;
    // printf("%e\n",data->pk[i]);
  }
}

void QuitCole05Pk(struct DataPk * data, struct DataPk * lin)
{
  free(data->k);
  free(data->pk);
  free(data->err);

  free(lin->k);
  free(lin->pk);
  free(lin->err);
}

double chi2_Cole05Pk(double *params, int Nparams, struct DataPk * data, struct DataPk * lin)
{
  // Make parameters readable
  double b0,A,Q;
  b0 = params[0];
  A = params[1];
  Q = params[2];

  // double b0,A,Q;
  // b0 = 1.0;
  // A = params[0];
  // Q = params[1];

  double Mod[data->Nk];
  double deviation[data->Nk];

  int i;
  double chi2 = 0.0;

  for (i=0; i<data->Nk; i++)
  {
    Mod[i] = b0*b0* ((1.0+Q*(data->k[i]*data->k[i])) / (1.0+A*data->k[i])) *lin->pk[i] ;
    deviation[i] = data->pk[i]-Mod[i];
    // printf("%e %e %e %e %e\n",lin->pk[i],Mod[i],data->pk[i],deviation[i],data->err[i]);
  }

  for(i=0; i<data->Nk; i++)
  {
    // for(j=0;j<N;j++)
    // {
      // chi2 += deviation[i] * D->InvC[i][j] * deviation[j];
      chi2 += deviation[i] * deviation[i] / (data->err[i] * data->err[i]) ;
    // }
  }

  if(!check_for_nan_and_inf(chi2,"Cole05pk chi2")) return -1;

  #ifdef DEBUG
  printf("chi2 = %e\n",chi2);
  #endif

  return chi2;
}
