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
#include <vectors.h>
#include <linalg.h>
#include <runtime_messages.h>
#include <file_manipulation.h>

double randrange(double minval, double maxval)
{
  return (minval + ((double)rand()/RAND_MAX)*(maxval-minval));
}

void updated_covariance(gsl_matrix *C, int Nrow, int Ncols, double **par ,int current_step,int Nparams)
{
  int i,j;
  double mean[Nparams]; for (i=0;i<Nparams; i++) mean[i] = 0.0;
  double val;

  for (j=0; j<current_step+1; j++)
  {
    for (i=0;i<Nparams; i++)
    {
      mean[i]+=par[j][i];
    }
  }
  for (i=0;i<Nparams; i++) mean[i] /= (double)(current_step+1);

  int k;
  for (i=0; i<Nparams; i++)
  {
    for (j=0;j<Nparams; j++)
    {
      val = 0.0;
      for (k=0; k<current_step+1;k++) val += (par[k][i]-mean[i])*(par[k][j]-mean[j]);
      gsl_matrix_set(C,i,j,val/((double)(current_step+1)));
    }
  }
}

double box_muller(double mu, double sigma)
{
  return (mu + sigma*( sqrt(-2.0*log(randrange(0.0,1.0))) * cos(2.0*M_PI*randrange(0.0,1.0)) ) );
}

void compute_proposed(double *theta_1,double*theta_0,double **par,
  gsl_matrix *current_covariance,gsl_matrix *step_cov,gsl_matrix *mixing,
  gsl_vector *diag_step_cov,int Nparams,int i_step)
{
  int i,j;
  double uncorrel[Nparams];

  if (i_step%1000==0.0)
  {
    // printf("ciao %i \n",i_step);
    updated_covariance(current_covariance,Nparams,Nparams,par,i_step,Nparams);

    gsl_matrix * cpy_step_cov = gsl_matrix_alloc(Nparams,Nparams);
    gsl_eigen_symmv_workspace * W = gsl_eigen_symmv_alloc (Nparams);

    // Compute step covariance matrix
    for (i=0; i<Nparams; i++)
    {
      for (j=0; j<Nparams; j++)
      {
        gsl_matrix_set(step_cov,i,j,NU*NU*gsl_matrix_get(current_covariance,i,j));
        gsl_matrix_set(cpy_step_cov,i,j,gsl_matrix_get(step_cov,i,j));
      }
    }

    // Compute mixing matrix
    // P | C = P L PT
    // L = digonal matrix given by the eigenvalues of C
    // P = matrix given by eigenvectors of C
    gsl_eigen_symmv (cpy_step_cov, diag_step_cov, mixing, W);

    // Compute uncorrelated gaussian sampling
    for (i=0; i<Nparams; i++)
    {
      uncorrel[i] = box_muller(0.0,sqrt(gsl_vector_get(diag_step_cov,i)));
    }

    gsl_matrix_free(cpy_step_cov);
    gsl_eigen_symmv_free(W);
  }
  else
  {
    // Compute uncorrelated gaussian sampling
    for (i=0; i<Nparams; i++)
    {
      uncorrel[i] = box_muller(0.0,sqrt(gsl_vector_get(diag_step_cov,i)));
    }
  }

  double incr[Nparams];
  for (i=0; i<Nparams; i++)
  {
    incr[i] = 0.0;
    for (j=0; j<Nparams; j++)
    {
      incr[i] += gsl_matrix_get(mixing,i,j) * uncorrel[j];
    }
  }

  for (i=0; i<Nparams; i++)
  {
    theta_1[i] = incr[i] + theta_0[i];
  }
}

void Gelman_Rubin(int Nchain, int Nparams, char chaindir[], char root_name[], char **paramnames)
{
  double param[Nparams];
  int indexchain,indexpar;

  int l,nl;
  int ChainLen[Nchain];
  double ParMeanTot[Nparams];
  double ParVarTot[Nparams];

  double **ParMean = allocate_double_matrix(Nchain,Nparams);
  double **ParVar = allocate_double_matrix(Nchain,Nparams);
  if(status!=SUCCESS) return;

  for (indexchain=0;indexchain<Nchain;indexchain++)
  {
    ChainLen[indexchain] = 0;
    for (indexpar=0;indexpar<Nparams;indexpar++)
    {
      ParMean[indexchain][indexpar] = 0.0;
      ParVar[indexchain][indexpar] = 0.0;
    }
  }
  for (indexpar=0;indexpar<Nparams;indexpar++)
  {
    ParMeanTot[indexpar] = 0.0;
    ParVarTot[indexpar] = 0.0;
  }

  double N;
  double chi2;
  FILE *f;

  char chainfile[STD_STR];
  for (indexchain=0;indexchain<Nchain;indexchain++)
  {
    sprintf(chainfile,"%s%s_%i.txt",chaindir,root_name,rank+1);
    nl = count_lines(chainfile);
    if(status!=SUCCESS) return;
    if(count_number_of_columns(chainfile,count_header_lines(chainfile))!=(2+Nparams))
    {
      frame("Are you sure your file has the correct layout?\n"
             "N -lnL/2 par1  ... parN\n"
             "Incorrect number of columns detected\n");
      status = IO_FAILURE;
      return;
    }
    f = fopen(chainfile,"r");
    for (l=0;l<nl;l++)
    {
      fscanf(f,"%lf %lf",&N,&chi2);
      for (indexpar=0;indexpar<Nparams;indexpar++)
      {
        fscanf(f,"%lf",&param[indexpar]);
        ParMean[indexchain][indexpar] += N*param[indexpar];
      }
      ChainLen[indexchain] += (int)N;
    }
    fclose(f);
  }

  int TotN = 0.0;
  for (indexchain=0;indexchain<Nchain;indexchain++) TotN += ChainLen[indexchain];

  for (indexpar=0;indexpar<Nparams;indexpar++)
  {
    for (indexchain=0;indexchain<Nchain;indexchain++)
    {
      ParMeanTot[indexpar] += ParMean[indexchain][indexpar];
    }
  }

  // Means
  for (indexchain=0;indexchain<Nchain;indexchain++)
  {
    for (indexpar=0;indexpar<Nparams;indexpar++)
    {
      ParMean[indexchain][indexpar] /= ChainLen[indexchain];
    }
  }
  for (indexpar=0;indexpar<Nparams;indexpar++)
  {
    ParMeanTot[indexpar]/=TotN;
  }


  for (indexchain=0;indexchain<Nchain;indexchain++)
  {
    sprintf(chainfile,"%s%s_%i.txt",chaindir,root_name,rank+1);
    nl = count_lines(chainfile);
    if(status!=SUCCESS) return;
    f = fopen(chainfile,"r");
    for (l=0;l<nl;l++)
    {
      fscanf(f,"%lf %lf",&N,&chi2);
      for (indexpar=0;indexpar<Nparams;indexpar++)
      {
        fscanf(f,"%lf",&param[indexpar]);
        ParVar[indexchain][indexpar] += N*pow(param[indexpar]-ParMean[indexchain][indexpar],2.0);
        ParVarTot[indexpar] += N*pow(param[indexpar]-ParMeanTot[indexpar],2.0);
      }
    }
    fclose(f);
  }

  for (indexchain=0;indexchain<Nchain;indexchain++)
  {
    for (indexpar=0;indexpar<Nparams;indexpar++)
    {
      ParVar[indexchain][indexpar]/=ChainLen[indexchain];
    }
  }
  for (indexpar=0;indexpar<Nparams;indexpar++)
  {
    ParVarTot[indexpar]/=TotN;
  }

  double W[Nparams];
  double B[Nparams];

  for (indexpar=0;indexpar<Nparams;indexpar++)
  {
    W[indexpar]=0.0;
    B[indexpar]=0.0;
    for (indexchain=0;indexchain<Nchain;indexchain++)
    {
      W[indexpar] += ParVar[indexchain][indexpar];
      B[indexpar] += (ChainLen[indexchain]/(Nchain-1.0))*pow(ParMean[indexchain][indexpar] - ParMeanTot[indexpar],2);
    }
    W[indexpar] /= Nchain;
  }

  double V[Nparams];
  double R[Nparams];
  for (indexpar=0;indexpar<Nparams;indexpar++)
  {
    V[indexpar] = ((ChainLen[0]-1.0)/ChainLen[0])*W[indexpar]
    + ((Nchain+1.0)/(Nchain*ChainLen[0]))*B[indexpar];
    R[indexpar] = V[indexpar]/W[indexpar];
  }

  printf("\nGelman-Rubin Convergence Test.\n\nR for each parameter is:\n");
  printf("\n  par\t\t\tR - 1\n\n");
  for (indexpar=0;indexpar<Nparams;indexpar++)
  {
    printf("  %s\t\t\t%lf\n",paramnames[indexpar],R[indexpar]-1.0);
  }
  printf("\n");

  deallocate_double_matrix(ParMean,Nchain,Nparams);
  deallocate_double_matrix(ParVar,Nchain,Nparams);
}
