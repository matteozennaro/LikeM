#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include <extern_global.h>
#include <file_manipulation.h>
#include <linalg.h>
#include <read_ini_file.h>
#include <runtime_messages.h>
#include <vectors.h>
#include <numerics.h>
#include <statistics.h>

#define LOW 1e-4
#define HIGH 1000
#define ABSPREC 0
#define RELPREC 1e-3

struct DataNonLocPk {
  int Nk;
  double *k;
  double *pk;
  double *err;
};

struct LinNonLocPk {
  int Nk;
  double *k;
  double *pk;
  double kmin;
  double kmax;
  double r;
  double know;
  gsl_interp_accel *acc;
  gsl_spline *splinePk;
};

void InitialiseNonLocRSPk(char parfile[], struct DataNonLocPk *gal, struct DataNonLocPk *data, struct LinNonLocPk *lin)
{
  char datagalfile[STD_STR];
  char dataccfile[STD_STR];
  char linfile[STD_STR];
  if(!read_string_from_file(parfile,"galpower",datagalfile)) read_err("galpower");
  if(!read_string_from_file(parfile,"nlpower",dataccfile)) read_err("nlpower");
  if(!read_string_from_file(parfile,"linpower",linfile)) read_err("linpower");
  if(!check_status(status)) mpi_exit(status);

  data->Nk = count_lines(dataccfile) - 1;
  data->k = allocate_double_vector(data->Nk);
  data->pk = allocate_double_vector(data->Nk);
  data->err = allocate_double_vector(data->Nk);
  int i;
  char buf[LONG_STR];
  double val,psn;
  FILE *f = fopen(dataccfile,"r");
  fscanf(f,"%lf %lf",&val,&psn);
  psn *= (2.0*M_PI)*(2.0*M_PI)*(2.0*M_PI);
  for(i=0;i<data->Nk;i++)
  {
    fscanf(f,"%lf %lf %lf %lf",&val,&data->k[i],&data->pk[i],&data->err[i]);
    data->pk[i] *= (2.0*M_PI)*(2.0*M_PI)*(2.0*M_PI);
    data->err[i] = sqrt(2.0/data->err[i])*data->pk[i];
    data->pk[i] -= psn;
  }
  fclose(f);

  gal->Nk = count_lines(datagalfile) - 1;
  gal->k = allocate_double_vector(gal->Nk);
  gal->pk = allocate_double_vector(gal->Nk);
  gal->err = allocate_double_vector(gal->Nk);
  double valgal,psngal;
  f = fopen(datagalfile,"r");
  fscanf(f,"%lf %lf",&valgal,&psngal);
  psn *= (2.0*M_PI)*(2.0*M_PI)*(2.0*M_PI);
  for(i=0;i<data->Nk;i++)
  {
    fscanf(f,"%lf %lf %lf %lf",&val,&gal->k[i],&gal->pk[i],&gal->err[i]);
    gal->pk[i] *= (2.0*M_PI)*(2.0*M_PI)*(2.0*M_PI);
    gal->err[i] = sqrt(2.0/gal->err[i])*gal->pk[i];
    gal->pk[i] -= psngal;
  }
  fclose(f);

  int hdlines = count_header_lines(linfile);
  lin->Nk = count_lines(linfile) - hdlines;
  lin->k = allocate_double_vector(lin->Nk);
  lin->pk = allocate_double_vector(lin->Nk);
  f = fopen(linfile,"r");
  for(i=0;i<hdlines;i++) fgets(buf,sizeof(buf),f);
  for(i=0;i<lin->Nk;i++)
  {
    fscanf(f,"%lf %lf",&lin->k[i],&lin->pk[i]);
  }
  fclose(f);

  lin->acc = gsl_interp_accel_alloc();
  lin->splinePk = gsl_spline_alloc(gsl_interp_cspline,lin->Nk);
  gsl_spline_init(lin->splinePk,lin->k,lin->pk,lin->Nk);
  lin->kmin = min_double_vec(lin->Nk,lin->k);
  lin->kmax = max_double_vec(lin->Nk,lin->k);
}

void QuitNonLocRSPk(struct DataNonLocPk *data, struct LinNonLocPk *lin)
{
  free(data->k);
  free(data->pk);
  free(data->err);

  free(lin->k);
  free(lin->pk);
  gsl_spline_free(lin->splinePk);
  gsl_interp_accel_free(lin->acc);
}

/******************************************************************************/

double P_at_k(double k, struct LinNonLocPk *p)
{
  double Pk = 0;
  if (k>p->kmin && k<p->kmax)
    Pk = gsl_spline_eval(p->splinePk,k,p->acc);
  else if(k<=p->kmin)
  {
    Pk = exp((log(p->pk[1]/p->pk[0])/log(p->k[1]/p->k[0]))*log(k/p->k[0]) + log(p->pk[0]));
  }
  else if(k>=p->kmax)
  {
    Pk = exp((log(p->pk[p->Nk-1]/p->pk[p->Nk-2])/log(p->k[p->Nk-1]/p->k[p->Nk-2]))*log(k/p->k[p->Nk-2]) + log(p->pk[p->Nk-2]));
  }
  else
  {
    printf("Unexpected error in interpolation!\n");
    status = NUM_FAILURE;
  }
  if(isnan(Pk))
  {
    printf("Pk is NaN in P_at_k!!! When this happened\nk = %e\nkmin = %e\nkmax = %e\n",
    k,p->kmin,p->kmax);
    status = NUM_FAILURE;
  }
  return Pk;
}

// generic power law interpolation of a power spectrum
double PSPL(double k, double *ko, double *pko, int N)
{
  double Pk = 0;
  double kmin = ko[0];
  double kmax = ko[N-1];

  if (k>kmin && k<kmax)
  {
    int i=0;
    while(ko[i] < k) i++;
    Pk = exp((log(pko[i]/pko[i-1])/log(ko[i]/ko[i-1]))*log(k/ko[i-1]) + log(pko[i-1]));
  }
  else if(k<=kmin)
  {
    Pk = exp((log(pko[1]/pko[0])/log(ko[1]/ko[0]))*log(k/ko[0]) + log(pko[0]));
  }
  else if(k>=kmax)
  {
    Pk = exp((log(pko[N-1]/pko[N-2])/log(ko[N-1]/ko[N-2]))*log(k/ko[N-2]) + log(pko[N-2]));
  }
  else
  {
    printf("Unexpected error in generic interpolation!\n");
    status = NUM_FAILURE;
  }
  if(isnan(Pk))
  {
    printf("Pk is NaN in P_at_k!!! When this happened\nk = %e\nkmin = %e\nkmax = %e\n",
    k,kmin,kmax);
    status = NUM_FAILURE;
  }
  return Pk;
}

/******************************************************************************/

double Pb2delta_x_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;
  double res;
  double r = p->r;
  double x = var;
  double k = p->know;

  double q = k*r;
  double sq_rt = sqrt(1.0+r*r-2.0*r*x);
  double modkq = k*sq_rt;
  double mu = (x-r)/sq_rt;
  double F2 = 5.0/7.0+0.5*mu*(1.0+2.0*r*(r-x))/(r*sq_rt) + 2.0*mu*mu/7.0;

  double Vol = (k*k*k*r*r/(4.0*M_PI*M_PI));

  res = Vol * P_at_k(q,p)*P_at_k(modkq,p) * F2;

  if(isnan(res))
  {
    printf("Pb2delta_x_integrand is returning a nan!\n");
    printf("k = %e\nr = %e\nx = %e\nmu = %e\nF2 = %e\nres = %e\n",k,r,x,mu,F2,res);
    status = NUM_FAILURE;
  }
  if(isinf(res))
  {
    printf("Pb2delta_x_integrand is returning a inf!\n");
    printf("k = %e\nr = %e\nx = %e\nmu = %e\nF2 = %e\nres = %e\n",k,r,x,mu,F2,res);
    status = NUM_FAILURE;
  }
  return res;
}

double Pb2delta_r_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;

  p->r = var;

  gsl_integration_cquad_workspace * wscq = gsl_integration_cquad_workspace_alloc(1000);

  double result, error;

  gsl_function F;
  F.function = &Pb2delta_x_integrand;
  F.params = p;

  size_t nevals;
  gsl_integration_cquad(&F, -1.0, 1.0, ABSPREC, RELPREC, wscq, &result, &error, &nevals);
  gsl_integration_cquad_workspace_free(wscq);

  if(isnan(result)) printf("Pb2delta_r_integrand is returning a nan!\n");
  if(isinf(result)) printf("Pb2delta_r_integrand is returning a inf!\n");

  return result;
}

/******************************************************************************/

double Pbs2delta_x_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;
  double res;
  double r = p->r;
  double x = var;
  double k = p->know;

  double q = k*r;
  double sq_rt = sqrt(1.0+r*r-2.0*r*x);
  double modkq = k*sq_rt;
  double mu = (x-r)/sq_rt;
  double F2 = 5.0/7.0+0.5*mu*(1.0+2.0*r*(r-x))/(r*sq_rt) + 2.0*mu*mu/7.0;
  double S2 = (x-r)*(x-r)/(1.0+r*r-2.0*r*x) - 1.0/3.0;

  double Vol = (k*k*k*r*r/(4.0*M_PI*M_PI));

  res = Vol * P_at_k(q,p)*P_at_k(modkq,p) * F2 * S2;

  if(isnan(res))
  {
    printf("Pbs2delta_x_integrand is returning a nan!\n");
    printf("k = %e\nr = %e\nx = %e\nmu = %e\nF2 = %e\nres = %e\n",k,r,x,mu,F2,res);
    status = NUM_FAILURE;
  }
  if(isinf(res))
  {
    printf("Pbs2delta_x_integrand is returning a inf!\n");
    printf("k = %e\nr = %e\nx = %e\nmu = %e\nF2 = %e\nres = %e\n",k,r,x,mu,F2,res);
    status = NUM_FAILURE;
  }
  return res;
}

double Pbs2delta_r_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;

  p->r = var;

  gsl_integration_cquad_workspace * wscq = gsl_integration_cquad_workspace_alloc(1000);

  double result, error;

  gsl_function F;
  F.function = &Pbs2delta_x_integrand;
  F.params = p;

  size_t nevals;
  gsl_integration_cquad(&F, -1.0, 1.0, ABSPREC, RELPREC, wscq, &result, &error, &nevals);
  gsl_integration_cquad_workspace_free(wscq);

  if(isnan(result)) printf("Pbs2delta_r_integrand is returning a nan!\n");
  if(isinf(result)) printf("Pbs2delta_r_integrand is returning a inf!\n");

  return result;
}

/******************************************************************************/

double Pb2s2delta_x_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;
  double res;
  double r = p->r;
  double x = var;
  double k = p->know;

  double q = k*r;
  double sq_rt = sqrt(1.0+r*r-2.0*r*x);
  double modkq = k*sq_rt;
  double S2 = (x-r)*(x-r)/(1.0+r*r-2.0*r*x) - 1.0/3.0;

  double Vol = (k*k*k*r*r/(4.0*M_PI*M_PI));

  res = -0.5* Vol * P_at_k(q,p) * ( 2.0/3.0*P_at_k(q,p) - P_at_k(modkq,p)*S2 );

  if(isnan(res))
  {
    printf("Pb2s2delta_x_integrand is returning a nan!\n");
    status = NUM_FAILURE;
  }
  if(isinf(res))
  {
    printf("Pb2s2delta_x_integrand is returning a inf!\n");
    status = NUM_FAILURE;
  }
  return res;
}

double Pb2s2delta_r_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;

  p->r = var;

  gsl_integration_cquad_workspace * wscq = gsl_integration_cquad_workspace_alloc(1000);

  double result, error;

  gsl_function F;
  F.function = &Pb2s2delta_x_integrand;
  F.params = p;

  size_t nevals;
  gsl_integration_cquad(&F, -1.0, 1.0, ABSPREC, RELPREC, wscq, &result, &error, &nevals);
  gsl_integration_cquad_workspace_free(wscq);

  if(isnan(result)) printf("Pb2s2delta_r_integrand is returning a nan!\n");
  if(isinf(result)) printf("Pb2s2delta_r_integrand is returning a inf!\n");

  return result;
}

/******************************************************************************/

double Pbs22delta_x_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;
  double res;
  double r = p->r;
  double x = var;
  double k = p->know;

  double q = k*r;
  double sq_rt = sqrt(1.0+r*r-2.0*r*x);
  double modkq = k*sq_rt;
  double S2 = (x-r)*(x-r)/(1.0+r*r-2.0*r*x) - 1.0/3.0;

  double Vol = (k*k*k*r*r/(4.0*M_PI*M_PI));

  res = -0.5* Vol * P_at_k(q,p) * ( 4.0/9.0*P_at_k(q,p) - P_at_k(modkq,p)*S2*S2 );

  if(isnan(res))
  {
    printf("Pbs22delta_x_integrand is returning a nan!\n");
    status = NUM_FAILURE;
  }
  if(isinf(res))
  {
    printf("Pbs22delta_x_integrand is returning a inf!\n");
    status = NUM_FAILURE;
  }
  return res;
}

double Pbs22delta_r_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;

  p->r = var;

  gsl_integration_cquad_workspace * wscq = gsl_integration_cquad_workspace_alloc(1000);

  double result, error;

  gsl_function F;
  F.function = &Pbs22delta_x_integrand;
  F.params = p;

  size_t nevals;
  gsl_integration_cquad(&F, -1.0, 1.0, ABSPREC, RELPREC, wscq, &result, &error, &nevals);
  gsl_integration_cquad_workspace_free(wscq);

  if(isnan(result)) printf("Pbs22delta_r_integrand is returning a nan!\n");
  if(isinf(result)) printf("Pbs22delta_r_integrand is returning a inf!\n");

  return result;
}

/******************************************************************************/

double Pb22delta_x_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;
  double res;
  double r = p->r;
  double x = var;
  double k = p->know;

  double q = k*r;
  double sq_rt = sqrt(1.0+r*r-2.0*r*x);
  double modkq = k*sq_rt;

  double Vol = (k*k*k*r*r/(4.0*M_PI*M_PI));

  res = -0.5* Vol * P_at_k(q,p) * ( P_at_k(q,p) - P_at_k(modkq,p) );

  if(isnan(res))
  {
    printf("Pb22delta_x_integrand is returning a nan!\n");
    status = NUM_FAILURE;
  }
  if(isinf(res))
  {
    printf("Pb22delta_x_integrand is returning a inf!\n");
    status = NUM_FAILURE;
  }
  return res;
}

double Pb22delta_r_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;

  p->r = var;

  gsl_integration_cquad_workspace * wscq = gsl_integration_cquad_workspace_alloc(1000);

  double result, error;

  gsl_function F;
  F.function = &Pb22delta_x_integrand;
  F.params = p;

  size_t nevals;
  gsl_integration_cquad(&F, -1.0, 1.0, ABSPREC, RELPREC, wscq, &result, &error, &nevals);
  gsl_integration_cquad_workspace_free(wscq);

  if(isnan(result)) printf("Pb22delta_r_integrand is returning a nan!\n");
  if(isinf(result)) printf("Pb22delta_r_integrand is returning a inf!\n");

  return result;
}

/******************************************************************************/

double Sigma_2_3_x_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;
  double res;
  double r = p->r;
  double x = var;
  double k = p->know;

  double q = k*r;
  double S2_1 = (x-r)*(x-r)/(1.0+r*r-2.0*r*x) - 1.0/3.0;
  double S2_2 = x*x - 1.0/3.0;

  double Vol = (k*k*k*r*r/(4.0*M_PI*M_PI));

  res = Vol * P_at_k(q,p) * ( 5.0/6.0 + 15.0/8.0*S2_1*S2_2 - 5.0/4.0*S2_1 );

  if(isnan(res))
  {
    printf("Sigma_2_3_x_integrand is returning a nan!\n");
    status = NUM_FAILURE;
  }
  if(isinf(res))
  {
    printf("Sigma_2_3_x_integrand is returning a inf!\n");
    status = NUM_FAILURE;
  }
  return res;
}

double Sigma_2_3_r_integrand(double var, void *params)
{
  struct LinNonLocPk *p = (struct LinNonLocPk *) params;

  p->r = var;

  gsl_integration_cquad_workspace * wscq = gsl_integration_cquad_workspace_alloc(1000);

  double result, error;

  gsl_function F;
  F.function = &Sigma_2_3_x_integrand;
  F.params = p;

  size_t nevals;
  gsl_integration_cquad(&F, -1.0, 1.0, ABSPREC, RELPREC, wscq, &result, &error, &nevals);
  gsl_integration_cquad_workspace_free(wscq);

  if(isnan(result)) printf("Sigma_2_3_r_integrand is returning a nan!\n");
  if(isinf(result)) printf("Sigma_2_3_r_integrand is returning a inf!\n");

  return result;
}

/******************************************************************************/

double compute_cquad(double know, double (* functionnow) (double var, void *params), struct LinNonLocPk *p)
{
  p->know = know;

  gsl_integration_cquad_workspace * wscq = gsl_integration_cquad_workspace_alloc(1000);

  double result, error;

  gsl_function F;
  F.function = functionnow;
  F.params = p;

  size_t nevals;
  gsl_integration_cquad(&F, LOW, HIGH, ABSPREC, RELPREC, wscq, &result, &error, &nevals);
  gsl_integration_cquad_workspace_free(wscq);

  if(isnan(result))
  {
    printf("Error! NaN found in cquad!\n");
    status = NUM_FAILURE;
  }

  return result;
}

double chi2_NonLocRSPk(double *params, int Nparams, struct DataNonLocPk * gal, struct DataNonLocPk * data, struct LinNonLocPk * lin)
{
  double b1, b2;
  double bs2, b3nl;
  b1 = params[0];
  b2 = params[1];

  bs2 = (-4.0/7.0)*(b1-1.0);
  b3nl = (32.0/315.0)*(b1-1.0);

  int i;

  int Nktab = 150;
  double *ktab = allocate_double_vector(Nktab);
  double *Pb2delta = allocate_double_vector(Nktab);
  double *Pbs2delta = allocate_double_vector(Nktab);
  double *Pb2s2delta = allocate_double_vector(Nktab);
  double *Pbs22delta = allocate_double_vector(Nktab);
  double *Pb22delta = allocate_double_vector(Nktab);
  double *Sigma_2_3 = allocate_double_vector(Nktab);

  for(i=0;i<Nktab;i++)
  {
    Pb2delta[i] = compute_cquad(ktab[i],&Pb2delta_r_integrand,lin);
    Pbs2delta[i] = compute_cquad(ktab[i],&Pbs2delta_r_integrand,lin);
    Pb2s2delta[i] = compute_cquad(ktab[i],&Pb2s2delta_r_integrand,lin);
    Pbs22delta[i] = compute_cquad(ktab[i],&Pbs22delta_r_integrand,lin);
    Pb22delta[i] = compute_cquad(ktab[i],&Pb22delta_r_integrand,lin);
    Sigma_2_3[i] = compute_cquad(ktab[i],&Sigma_2_3_r_integrand,lin);
  }

  if(status!=SUCCESS) return -1;

  double Mod;
  double diff;

  double chi2 = 0.0;

  for(i=0; i<data->Nk; i++)
  {
    Mod = b1*b1 * data->pk[i]
      + 2.0*b1*b2*PSPL(data->k[i],ktab,Pb2delta,Nktab)
      + 2.0*bs2*b1*PSPL(data->k[i],ktab,Pbs2delta,Nktab)
      + b2*b2*PSPL(data->k[i],ktab,Pb22delta,Nktab)
      + 2.0*b2*bs2*PSPL(data->k[i],ktab,Pb2s2delta,Nktab)
      + bs2*bs2*PSPL(data->k[i],ktab,Pbs22delta,Nktab)
      + 2.0*b1*b3nl*PSPL(data->k[i],ktab,Sigma_2_3,Nktab)*P_at_k(data->k[i],lin);

    diff = (gal->pk[i]-Mod) / gal->err[i];
    chi2 += diff*diff;
  }

  if(!check_for_nan_and_inf(chi2,"Cole05pk chi2")) return -1;

  #ifdef DEBUG
  printf("chi2 = %e\n",chi2);
  #endif

  free(ktab);
  free(Pb2delta);
  free(Pbs2delta);
  free(Pb2s2delta);
  free(Pbs22delta);
  free(Pb22delta);
  free(Sigma_2_3);

  return chi2;
}
