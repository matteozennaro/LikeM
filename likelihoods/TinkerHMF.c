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

struct HMFDATA {
	int N;
	double * BinMin;
	double * BinMax;
	double * BinMean;
	double * HMF;
	double * HMFErr;
	double **C;
	double **InvC;
	double Mmin;
	double Mmax;
	double hSim;
	double VolSim;
	double Omega_m;
	double rho_bg;
};

struct STab {
	int N;
	double Mmin;
  double Mmax;
  double * M;
  double * Sigma;
};

void read_histo(char file[],int N,double *min,double *max,
	double *mean,double *mf,double *err)
{
	int num_hd = count_header_lines(file);
	if(status!=SUCCESS) return;
	int num_col = count_number_of_columns(file,num_hd);
	if(status!=SUCCESS) return;

	if (num_col!=5)
	{
		char error[LONG_STR];
		sprintf(error,"Error! Unexpected number of columns in file\n%s\n",file);
		frame(error);
		status = IO_FAILURE;
		return;
	}

  FILE *f = fopen(file,"r");
  if (f==NULL)
  {
    printf("\n\tError opening file %s\n",file);
    status = IO_FAILURE;
  }
  else
  {
    int i;
		char buf[LONG_STR];
		for (i=0;i<num_hd;i++) if(fgets(buf,sizeof(buf),f)==NULL) exit(-1);
    for (i=0;i<N;i++)
		{
			fscanf(f,"%lf %lf %lf %lf %lf",
			&min[i],&max[i],&mean[i],&mf[i],&err[i]);
		}
		fclose(f);
  }
}


void read_hmf_cov(char covf[], int N, double **C)
{
	FILE *f = fopen(covf,"r");
	if(f==NULL)
	{
		printf("\n\tError opening file %s\n",covf);
    status = IO_FAILURE;
	}
	else
	{
		int i,j;
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				fscanf(f,"%lf",&C[i][j]);
			}
		}
		fclose(f);
	}
}

void compute_inverse_covariance(double **C, double **InvC, int N)
{
	int i,j;
	gsl_matrix * cov = gsl_matrix_alloc(N,N);
	gsl_matrix * inv = gsl_matrix_alloc(N,N);

	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			gsl_matrix_set(cov,i,j,C[i][j]);
		}
	}

	// invert covariance matrix
	int signum;
	gsl_permutation * perm = gsl_permutation_alloc(N);
	gsl_linalg_LU_decomp(cov,perm,&signum);
	gsl_linalg_LU_invert(cov,perm,inv);

	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			InvC[i][j] = gsl_matrix_get(inv,i,j);
		}
	}

	#ifdef DEBUG
	if(rank==0)
	{
		double **I = allocate_double_matrix(N,N);
		if(status!=SUCCESS) return;
		multiply_squared_matrices(C,InvC,I,N);
		printf("Check that covariance matrix inversion was successful:\n\n");
		for(i=0;i<N;i++)
		{
			for(j=0;j<N;j++)
			{
				printf("%e ",I[i][j]);
			}
			printf("\n");
		}
		printf("\n");

		deallocate_double_matrix(I,N,N);
	}
	#endif

	gsl_matrix_free(cov);
	gsl_matrix_free(inv);
	gsl_permutation_free(perm);
}

void FillSigmaTab(char file[], struct STab * SigmaTab)
{
	int num_hd = count_header_lines(file);
	if(status!=SUCCESS) return;
	int num_col = count_number_of_columns(file,num_hd);
	if(status!=SUCCESS) return;

	if (num_col!=2)
	{
		char error[LONG_STR];
		sprintf(error,"Error! Unexpected number of columns in file\n%s\n",file);
		frame(error);
		status=IO_FAILURE;
		return;
	}

  FILE *f = fopen(file,"r");
  if (f==NULL)
  {
    printf("Error, file %s not found.\n",file);
    status=IO_FAILURE;
  }
  else
  {
    int i;
		char buf[2000];
		for (i=0;i<num_hd;i++) if(fgets(buf,sizeof(buf),f)==NULL) exit(-1);
    for (i=0;i<SigmaTab->N;i++)
    {
      fscanf(f,"%lf %lf",&SigmaTab->M[i],&SigmaTab->Sigma[i]);
    }
    fclose(f);
  }
}

void FreeSigmaTab(struct STab * SigmaTab)
{
  free(SigmaTab->M);
  free(SigmaTab->Sigma);
}

void FreeHMFdata(struct HMFDATA * D)
{
  free(D->BinMin);
  free(D->BinMax);
  free(D->BinMean);
  free(D->HMF);
	free(D->HMFErr);
	deallocate_double_matrix(D->C,D->N,D->N);
	deallocate_double_matrix(D->InvC,D->N,D->N);
}

double pkinterp(int n, double *x, double *y, double x0)
{
  int i=0;

  if (x0 > x[n-1])
  {
    i=n-1;
    return ((y[i]-y[i-1])/(x[i]-x[i-1])*x0 + y[i] - (y[i]-y[i-1])/(x[i]-x[i-1])*x[i]);
  }
  else if (x0 < x[0])
  {
    i=1;
    return ((y[i]-y[i-1])/(x[i]-x[i-1])*x0 + y[i] - (y[i]-y[i-1])/(x[i]-x[i-1])*x[i]);
  }
  else
  {
    while (x0 > x[i]) i++;
    if (i==0) return y[i];
    else return ((y[i]-y[i-1])/(x[i]-x[i-1])*x0 + y[i] - (y[i]-y[i-1])/(x[i]-x[i-1])*x[i]);
  }
}

double W2tophat(double k, double R)
{
  double W = 3.0*(sin(k*R)-(k*R)*cos(k*R))/((k*R)*(k*R)*(k*R));
  return W*W;
}

double func(double lnk, double lnP, double R)
{
  double k=exp(lnk);
  double pk=exp(lnP);
  return ((k*k*k*pk*W2tophat(k,R))/(2.*M_PI*M_PI));
  // return ((k*k*k*pk*W2(k,R)));
}

double Sigma(double M, int k_num, double *k, double *Pk, double rho_bg)
{
  double R = pow((3.0*M)/(4.0*rho_bg*M_PI),1.0/3.0);

  double kmin, kmax, kint;

  kmin = min_double_vec(k_num,k);
  kmax = max_double_vec(k_num,k);

  int n_step = 100000;

  double step = (kmax-kmin)/(n_step-1);

  double Integral = 0.0;
  double K0,K1,P0,P1;

  for (kint=kmin; kint<kmax; kint+=step)
  {
    K0 = kint;
    P0 = pkinterp(k_num,k,Pk,K0);
    K1 = kint + step;
    P1 = pkinterp(k_num,k,Pk,K1);
    Integral += 0.5*(func(K0,P0,R)+func(K1,P1,R))*step;
  }
  return sqrt(Integral);
}

double SigmaInterp(double M, struct STab * S)
{
  int i=0;
  if (M > S->Mmax)
  {
    frame("Error!! Requested sigma for a M above the maximum\n"
    "mass in pre-computed table.\n");
		status = FAILURE;
		return -1;
  }
  else if (M < S->Mmin)
  {
    frame("Error!! Requested sigma for a M below the minimum\n"
    "mass in pre-computed table.\n");
		status = FAILURE;
		return -1;
  }
  else
  {
    while(M > S->M[i]) i++;
    double m = (S->Sigma[i]-S->Sigma[i-1])/(S->M[i]-S->M[i-1]);
    return S->Sigma[i-1] + m*(M - S->M[i-1]);
  }
}

double ComputeNu(double M, struct STab * S)
{
	double sM = SigmaInterp(M,S);
	if(status!=SUCCESS) return -1;
	double dc = 1.686;
	return (dc*dc)/(sM*sM);
}

void ComputeSigmaTab(char pkf[], char sigmaf[], int N, double Mmin, double Mmax, double rho_bg)
{
  int i;
  double lMMin = log10(Mmin);
  double lMMax = log10(Mmax);
  double lMStep = (lMMax-lMMin)/(N-1);

	int headlines = count_header_lines(pkf);
	if(status!=SUCCESS) return;
	int k_num = count_lines(pkf) - headlines;
	if(status!=SUCCESS) return;
	double *k_in = allocate_double_vector(k_num);
	if(status!=SUCCESS) return;
	double *Pk_in = allocate_double_vector(k_num);
	if(status!=SUCCESS) return;

	FILE *f;
	f = fopen(pkf,"r");
	if(f!=NULL)
	{
		char buf[LONG_STR];
		for(i=0;i<headlines;i++) fgets(buf,sizeof(buf),f);
		for(i=0;i<k_num;i++)
		{
			fscanf(f,"%lf %lf",&k_in[i],&Pk_in[i]);
			// k_in[i] = log(k_in[i]);
			// Pk_in[i] = log(Pk_in[i]);
		}
		fclose(f);
	}
	else
	{
		char message[LONG_STR];
		sprintf(message,"Error! Unable to open file %s\n",pkf);
		frame(message);
		status = IO_FAILURE;
		return;
	}

  time_t start=time(NULL);
  printf("Computing sigma table...\n");
	printf("File : %s\n",pkf);
	printf("Header lines : %i\n",headlines);
	printf("Begins with: %e %e\n",k_in[0],Pk_in[0]);

  f = fopen(sigmaf,"w");
	if(f!=NULL)
	{
		double M,R,Sinterp;
	  for (i=0;i<N;i++)
	  {
			M = pow(10.0,lMMin+i*lMStep);
			R = pow((3.0*M)/(4.0*rho_bg*M_PI),1.0/3.0);
			Sinterp = sqrt(SIGMA2(k_num,k_in,Pk_in,R));
			if(!check_for_nan_and_inf(Sinterp,"interpolated sigma")) return;
			// fprintf(f,"%.10e\t%.10e\n",M,Sigma(M,k_num,k_in,Pk_in,rho_bg));
	    fprintf(f,"%.10e\t%.10e\n",M,Sinterp);
	  }
  	fclose(f);
	}
	else
	{
		char message[LONG_STR];
		sprintf(message,"Error! Unable to create file %s\n",sigmaf);
		frame(message);
		status = IO_FAILURE;
	}
  time_t stop = time(NULL);
  printf("Done, %ld seconds elapsed.\n",stop-start);
	free(k_in);
	free(Pk_in);
}

void InitialiseTinker(char ini[], struct HMFDATA * D, struct STab * S)
{
	char histof[LONG_STR];
  char covf[LONG_STR];
  char compute_s;
  char pkf[LONG_STR];
  char sigmaf[LONG_STR];

	if(!read_string_from_file(ini,"histo_file",histof))
	{read_err("histo_file"); status=IO_FAILURE; return;}
	if(!read_string_from_file(ini,"covariance",covf))
	{read_err("covariance"); status=IO_FAILURE; return;}
	if(!read_string_from_file(ini,"pk_linear",pkf))
	{read_err("pk_linear");status=IO_FAILURE; return;}
	if(!read_string_from_file(ini,"sigma_file",sigmaf))
	{read_err("sigma_file");status=IO_FAILURE; return;}
	if(!read_bool_from_file(ini,"compute_sigma_tab",&compute_s))
	{read_err("compute_sigma_tab");status=IO_FAILURE; return;}

  D->N = count_lines(histof)-count_header_lines(histof);
	if(status!=SUCCESS) return;
  D->BinMin = allocate_double_vector(D->N);
  D->BinMax = allocate_double_vector(D->N);
  D->BinMean = allocate_double_vector(D->N);
  D->HMF = allocate_double_vector(D->N);
	D->HMFErr = allocate_double_vector(D->N);
	D->C = allocate_double_matrix(D->N,D->N);
	D->InvC = allocate_double_matrix(D->N,D->N);
	if(status!=SUCCESS) return;
  D->hSim = 0.67;
	D->Omega_m = 0.32;
  D->VolSim = pow(2000.0,3);
	D->rho_bg = D->Omega_m * (3.0*100.0*100.0)/(8.0*M_PI*4.301e-9);
  read_histo(histof,D->N,D->BinMin,D->BinMax,D->BinMean,D->HMF,D->HMFErr);
	if(status!=SUCCESS) return;
	read_hmf_cov(covf,D->N,D->C);
	if(status!=SUCCESS) return;
	compute_inverse_covariance(D->C,D->InvC,D->N);
	if(status!=SUCCESS) return;
  D->Mmin = min_double_vec(D->N,D->BinMin);
  D->Mmax = max_double_vec(D->N,D->BinMax);

	S->Mmin = D->Mmin - 0.1*D->Mmin;
  S->Mmax = D->Mmax + 0.1*D->Mmax;

	if(compute_s=='T' && rank==0)
	{
		ComputeSigmaTab(pkf,sigmaf,300,S->Mmin,S->Mmax,D->rho_bg);
	}
	if(!check_status(status)) return;
	MPI_Barrier(MPI_COMM_WORLD);

	S->N = count_lines(sigmaf)-count_header_lines(sigmaf);
	if(status!=SUCCESS) return;
  S->M = allocate_double_vector(S->N);
  S->Sigma = allocate_double_vector(S->N);
	if(status!=SUCCESS) return;
  FillSigmaTab(sigmaf,S);
	if(status!=SUCCESS) return;
}

void QuitTinker(struct HMFDATA * HMFdata,struct STab * SigmaTab)
{
  FreeHMFdata(HMFdata);
  FreeSigmaTab(SigmaTab);
}

double chi2_Tinker(double *params, int Nparams, struct HMFDATA * D,
	struct STab * S)
{
  // Make parameters readable
  double A,a,b,c; //d;
  A = params[0];
  a = params[1];
  b = params[2];
  c = params[3];
	// d = params[4];
	// double G = 4.301e-9; //km2 Mpc MSun-1 s-2
  // double rho_c = (3.0*100.0*100.0)/(8.0*M_PI*G);
  // double rho_bg = 0.32*rho_c;
	double rho_bg = D->rho_bg;

  int i,j;
  int N = D->N;
	double sm,Ml,Mr,Deriv,f;
	double TOL = 0.01;
  double Mod[N];
  double deviation[N];

  double chi2 = 0.0;

  for (i=0;i<N;i++)
  {
		sm = SigmaInterp(D->BinMean[i],S);
		if(status!=SUCCESS) return -1;
		Ml = D->BinMean[i]*(1.0-TOL);
		Mr = D->BinMean[i]*(1.0+TOL);
		Deriv = log((1./SigmaInterp(Mr,S)) / (1./SigmaInterp(Ml,S))) / (Mr - Ml);
		if(status!=SUCCESS) return -1;
		f = A * (pow(sm/b,-a) + 1.0) * exp(-c/(sm*sm));
		Mod[i] = f * (rho_bg/D->BinMean[i]) * Deriv;

    deviation[i] = D->HMF[i]-Mod[i];
    // chi2 += ((deviation*deviation)/(HMFdata->HMFErr[i]*HMFdata->HMFErr[i]));
		// printf("%e %e %e %e\n",HMFdata->BinMean[i],HMFdata->HMF[i],HMFdata->HMFErr[i],Mod);
  }

	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			chi2 += deviation[i] * D->InvC[i][j] * deviation[j];
		}
	}

	// printf("%i, %e\n",rank,chi2);

	if(!check_for_nan_and_inf(chi2,"Tinker chi2")) return -1;

  #ifdef DEBUG
  printf("chi2 = %e\n",chi2);
  #endif

  return chi2;
}
