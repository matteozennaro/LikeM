#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

extern int error_mode;
extern char classdir[500];
extern char dir[500];

void set_planck_bf(double *planck_bf, char bestfit_st[])
{
  if (strcmp(bestfit_st,"planck")==0)
  {
    // Planck Best Fit
    planck_bf[0] = 0.02225 ; //Ob0h2
    planck_bf[1] = 0.1198  ; //Oc0h2
    planck_bf[2] = 0.079   ; //tau
    planck_bf[3] = 0.06    ; //mnu
    planck_bf[4] = 3.094   ; //log10^10As
    planck_bf[5] = 0.9645  ; //ns
    planck_bf[6] = 67.27   ; // H0
  }
  else if (strcmp(bestfit_st,"sim")==0)
  {
    // Sim Best Fit
    planck_bf[0] = 0.022445000000000007 ; //Ob0h2
    // planck_bf[1] = 0.11798204230191113  ; //Oc0h2
    planck_bf[1] = 0.11937779063774966  ; //Oc0h2
    planck_bf[2] = 0.0925 ;               //tau
    // planck_bf[3] = 0.30   ;               //mnu
    planck_bf[3] = 0.17  ;
    planck_bf[4] = 3.0570625287055968   ; //log10^10As
    planck_bf[5] = 0.96  ;                //ns
    planck_bf[6] = 67.0  ;                // H0
  }
  else
  {
    printf("Error: No legal choice for bestfit\n");
    error_mode=-1;
    return;
  }
}

double chi2_planck(double * theta_orig, double * theta_bf,
                   gsl_matrix *invC, int Nparams)
{
  double chi2=0.0;

  // REORDER TO AGREE WITH PLANCK ORDER
  double theta[Nparams];
  // h,OB0h2,OC0h2,mnu,tau,As,ns
  // 0 1    2      3   4   5  6

  theta[0] = theta_orig[1];
  theta[1] = theta_orig[2];
  theta[2] = theta_orig[4];
  theta[3] = theta_orig[3];
  theta[4] = log(10.0*theta_orig[5]);
  theta[5] = theta_orig[6];
  theta[6] = 100.0*theta_orig[0];

  int i,j;
  for (i=0; i<Nparams; i++)
  {
    for (j=0; j<Nparams; j++)
    {
      chi2 += (theta[i] - theta_bf[i]) * gsl_matrix_get(invC,i,j) * (theta[j] - theta_bf[j]);
    }
  }

  #ifdef DEBUG
    for (i=0; i<Nparams; i++)
    {
      for (j=0; j<Nparams; j++)
      {
        printf("%.2e  ",gsl_matrix_get(invC,i,j));
      }
      printf("\n");
    }
    printf("\n");
    printf("chi2 planck: %.10e\n",chi2);
  #endif

  return chi2;
}
