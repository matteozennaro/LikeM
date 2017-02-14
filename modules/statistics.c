#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <runtime_messages.h>
#include <vectors.h>
#include <file_manipulation.h>
#include <extern_global.h>
#include <numerics.h>

double W2(double k, double R)
{
  double W = 3.0*(sin(k*R)-(k*R)*cos(k*R))/((k*R)*(k*R)*(k*R));
  return W*W;
}

double func1(double lnk, double lnP, double r, double R)
{
  double k=exp(lnk);
  double pk=exp(lnP);
  return ((k*k*k*pk*W2(k,R)*(sin(k*r)/(k*r)))/(2.0*M_PI*M_PI));
}

double func2(double lnk, double lnP, double R)
{
  double k=exp(lnk);
  double pk=exp(lnP);
  return ((k*k*k*pk*W2(k,R))/(2.0*M_PI*M_PI));
}

double ETA(int k_num, double *k_o, double *Pk_o, double r, double R)
{
  int i;

  double *k = allocate_double_vector(k_num);
  double *Pk = allocate_double_vector(k_num);
  if(status!=SUCCESS) return -1;

  for (i=0; i<k_num; i++) k[i] = log(k_o[i]);
  for (i=0; i<k_num; i++) Pk[i] = log(Pk_o[i]);

  double kmin, kmax, kint;

  kmin = min_double_vec(k_num,k);
  kmax = max_double_vec(k_num,k);

  int n_step = 100000;

  double step = (kmax-kmin)/(n_step-1);

  double Integral1 = 0.0;
  double Integral2 = 0.0;
  double K0,K1,P0,P1;

  for (kint=kmin; kint<kmax; kint+=step)
  {
    K0 = kint;
    P0 = lin_interp(K0,k,Pk,k_num);
    K1 = kint + step;
    P1 = lin_interp(K1,k,Pk,k_num);

    Integral1 += 0.5*(func1(K0,P0,r,R)+func1(K1,P1,r,R))*step;
    Integral2 += 0.5*(func2(K0,P0,R)+func2(K1,P1,R))*step;
  }
  #ifdef DEBUG
  printf("\nR in func = %lf\n",R);
  printf("Eta in func = %e\n",Integral1/Integral2);
  #endif

  deallocate_double_vector(k);
  deallocate_double_vector(Pk);

  return Integral1/Integral2;
}

double SIGMA2(int k_num, double *k_o, double *Pk_o, double R)
{
  int i;

  double *k = allocate_double_vector(k_num);
  double *Pk = allocate_double_vector(k_num);
  if(status!=SUCCESS) return -1;

  for (i=0; i<k_num; i++) k[i] = log(k_o[i]);
  for (i=0; i<k_num; i++) Pk[i] = log(Pk_o[i]);

  for (i=0; i<k_num; i++)
  {
    if(!check_for_nan_and_inf(k[i],"k value")) return -1;
    if(!check_for_nan_and_inf(Pk[i],"Pk value")) return -1;
  }

  double kmin, kmax, kint;

  kmin = min_double_vec(k_num,k);
  kmax = max_double_vec(k_num,k);

  //printf("\tk_min = %.8lf\n",exp(kmin));
  //printf("\tk_max = %.8lf\n",exp(kmax));

  int n_step = 100000;

  double step = (kmax-kmin)/(n_step-1);

  double Integral2 = 0.0;
  double K0,K1,P0,P1;

  for (kint=kmin; kint<kmax; kint+=step)
  {
    K0 = kint;
    P0 = lin_interp(K0,k,Pk,k_num);
    K1 = kint + step;
    P1 = lin_interp(K1,k,Pk,k_num);

    Integral2 += 0.5*(func2(K0,P0,R)+func2(K1,P1,R))*step;
  }

  deallocate_double_vector(k);
  deallocate_double_vector(Pk);

  return Integral2;
}

double XI(int k_num, double *k_o, double *Pk_o, double r, double R)
{
  int i;

  double *k = allocate_double_vector(k_num);
  double *Pk = allocate_double_vector(k_num);
  if(status!=SUCCESS) return -1;

  for (i=0; i<k_num; i++) k[i] = log(k_o[i]);
  for (i=0; i<k_num; i++) Pk[i] = log(Pk_o[i]);

  double kmin, kmax, kint;

  kmin = min_double_vec(k_num,k);
  kmax = max_double_vec(k_num,k);

  int n_step = 100000;

  double step = (kmax-kmin)/(n_step-1);

  double Integral1 = 0.0;
  double K0,K1,P0,P1;

  for (kint=kmin; kint<kmax; kint+=step)
  {
    K0 = kint;
    P0 = lin_interp(K0,k,Pk,k_num);
    K1 = kint + step;
    P1 = lin_interp(K1,k,Pk,k_num);

    Integral1 += 0.5*(func1(K0,P0,r,R)+func1(K1,P1,r,R))*step;
  }

  deallocate_double_vector(k);
  deallocate_double_vector(Pk);

  return Integral1;
}
