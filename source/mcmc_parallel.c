#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <mpi.h>

#include <define_global.h>

struct HMFDATA {
	int N;
	double * BinMin;
	double * BinMax;
	double * BinMean;
	double * HMF;
	double * HMFErr;
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

#include <general_purpose.h>
#include <stat.h>
#include <read_ini_file.h>
#include <TinkerHMF.h>

int main (int argc, char * argv[])
{
	int ierr, size, rank;

	ierr = MPI_Init(&argc,&argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank==0)
	{
		// STARTING EVERYTHING...
		if (argc == 2)
		{
			if (strcmp(argv[1],"--help")==0)
			{
				print_help();
				exit(-1);
			}
			else
			{
				printf("\n*************************************************************\n");
				printf("\n                Welcome in mcmc!\n\n");
				printf("*************************************************************\n");

				getcwd(dir,sizeof(dir));
				srand(time(NULL));
				#ifdef DEBUG
					srand(5);
				#endif
				verb=1;

				int proc;
				int seed;
				for(proc=1; proc<size; proc++)
				{
					seed = rand();
					MPI_Send(&seed,1,MPI_INT,proc,0,MPI_COMM_WORLD);
				}
			}
		}
		else
		{
			printf("\n\tError! You must specify an input parameter file\n");
			printf("\tFor further info on the usage of mcmc, please type\n");
			printf("\t\t./mcmc --help\n\n");
			error_mode = -1;
			return error_mode;
		}
	}
	else
	{
		int proc;
		int seed;
		for(proc=1; proc<size; proc++)
		{
			if(rank==proc) MPI_Recv(&seed,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		srand(seed);
	}


	read_param_file(argv[1],rank);
	printf("\n");

	NU = 2.38/sqrt(Nparams);

	// PREPARE HALO MASS FUNCTION
	struct HMFDATA HMFdata;
	struct STab SigmaTab;
	InitialiseTinker(histofile,mfcovfile,&HMFdata,compute_sigma_tab,pk_file,sigmafile,&SigmaTab);

	int Nchain = size;
	int chain_ind = rank;

	/****************************************************************************/
	// DECLARATIONS
	int i_step,i_par;
	int N_accepted,N_tot;
	double param[Nparams];
	double par_step[Nparams];
	double theta_0[Nparams];
	double theta_1[Nparams];
	double chi2_0,chi2_1;
	double acceptance;

	int Nrow,Ncols;
	gsl_matrix *current_covariance;

	double **par;

	char outfile[500];
	time_t start,stop;

	for(i_par=0;i_par<Nparams;i_par++)
	{
		param[i_par] = param_read[i_par];
		par_step[i_par] = par_step_read[i_par];
	}

	/****************************************************************************/

	#ifdef DEBUG
		Nstep = 3;
		int i;
		for (i=0; i<Nparams; i++)
		{
			printf("%.10e	%.10e\n",param[i],par_step[i]);
		}
	#endif

	/****************************************************************************/
	// RANDOMIZE STARTING POINT
	for (i_par=0; i_par<Nparams; i_par++)
	{
		redoparam:
		param[i_par] = randrange(param[i_par]-5.0*par_step[i_par],
														 param[i_par]+5.0*par_step[i_par]);
		if (param[i_par] < 0.0) goto redoparam;
	}

	#ifdef DEBUG
		Nstep = 3;
		for (i=0; i<Nparams; i++)
		{
			printf("%.10e	%.10e\n",param[i],par_step[i]);
		}
	#endif

	/****************************************************************************/
	// METROPOLIS - HASTINGS
	start = time(NULL);

	N_accepted = N_tot = 0;
	sprintf(outfile,"%s_%i_raw.txt",root_name,chain_ind+name_offset); // printf("%s\n\n",outfile);

	i_step=0;
	printf("Starting chain %i!\n",chain_ind+name_offset);

	for (i_par=0; i_par<Nparams; i_par++) theta_0[i_par] = param[i_par];

	// redo1:
	chi2_0 = chi2_Tinker(theta_0,Nparams,&HMFdata,&SigmaTab);
	if (error_mode==-1) {printf("Internal error!\n"); goto wayout;}

	if (append == 1)
	{
		// retrieve last point from previously computed chain
		int nl = count_lines(outfile);
		double val;
		FILE * last = fopen(outfile,"r");
		for (i_step=0; i_step<nl; i_step++)
		{
			fscanf(last,"%lf",&val);
		for (i_par=0; i_par<Nparams; i_par++)
		{
			fscanf(last,"%lf",&theta_0[i_par]);
		}
		}
		fclose(last);

		printf("Last line in chain %i\n",chain_ind+name_offset);
		for (i_par=0; i_par<Nparams; i_par++) printf("%e  ",theta_0[i_par]);
		printf("\n");
	}

	Nrow=Ncols=Nparams;
	current_covariance = gsl_matrix_alloc(Nrow,Ncols);
	par = allocate_double_matrix(Nstep,Nparams);

	FILE * out;
	if (append==1) out=fopen(outfile,"a");
	else out=fopen(outfile,"w");

	for (i_par=0; i_par<Nparams; i_par++) par[0][i_par] = theta_0[i_par];
	for (i_par=0; i_par<Nparams; i_par++) gsl_matrix_set(current_covariance,i_par,i_par,par_step[i_par]*par_step[i_par]);

	int ii,jj;

	gsl_matrix * mixing = gsl_matrix_alloc(Nparams,Nparams);
	gsl_matrix * step_cov = gsl_matrix_alloc(Nparams,Nparams);
	gsl_vector * diag_step_cov = gsl_vector_alloc(Nparams);

	for (i_step=0; i_step<Nstep; i_step++)
	{
		// redo2:
		if (i_step>999)
		{
			compute_proposed(theta_1,theta_0,par,current_covariance,step_cov,mixing,diag_step_cov,Nparams,i_step);
		}
		else
		{
			for (i_par=0; i_par<Nparams; i_par++) theta_1[i_par] = theta_0[i_par] + par_step[i_par]*randrange(-1.0,1.0);
		}
		// for (i_par=0; i_par<Nparams; i_par++)
		//  if(theta_1[i_par]<0.0) goto redo2;


		chi2_1 = chi2_Tinker(theta_1,Nparams,&HMFdata,&SigmaTab);
		if (error_mode==-1) {printf("Internal error!\n"); goto wayout;}

		acceptance = chi2_1 - chi2_0;

		if (acceptance < -2.0*log(randrange(0.0,1.0)))
		{
			for (i_par=0; i_par<Nparams; i_par++) theta_0[i_par] = theta_1[i_par];
			chi2_0 = chi2_1;
			N_accepted++;
		}

		N_tot++;

		fprintf(out,"%.10e\t",chi2_0);
		for (i_par=0; i_par<Nparams; i_par++) fprintf(out,"%.10e\t",theta_0[i_par]);
		fprintf(out,"\n");
		fflush(out);

		for (i_par=0; i_par<Nparams; i_par++) par[i_step][i_par] = theta_0[i_par];
	}
	fclose(out);
	stop = time(NULL);
	printf("\n");

	printf("\nCovariance from chain %i:\n\n",chain_ind);
	for (ii=0;ii<Nparams;ii++)
	{
		for(jj=0;jj<Nparams;jj++) printf("%.3e  ",gsl_matrix_get(current_covariance,ii,jj));
		printf("\n");
	}
	printf("\n");

	printf("Ch %i, Number of accepted steps = %i\n"
				 "Ch %i, Acceptance rate = %lf\n"
				 "Ch %i, Total elapsed time = %ld s.\n",
				 chain_ind+name_offset,N_accepted,
				 chain_ind+name_offset,(double)N_accepted/N_tot,
				 chain_ind+name_offset,stop-start);
	/****************************************************************************/
	// REARRANGE OUTPUT FILE FOR USING getdist
	output_for_getdist(chain_ind+name_offset);

	// IF REQUESTED, PLOT CHAIN
	// python_plot_chain(root_name);

	/****************************************************************************/
	// PREPARING TO QUIT
	wayout:
	QuitTinker(&HMFdata,&SigmaTab);
	gsl_matrix_free(current_covariance);
	gsl_matrix_free(mixing);
	gsl_matrix_free(step_cov);
	gsl_vector_free(diag_step_cov);
	// gsl_permutation_free(PERM);
	free_double_matrix(par,Nstep,Nparams);

	// free_double_vector(ytab);
	// free_double_vector(FFtab);
	// free_double_vector(GGtab);

	// COMPUTE CONVERGENCE STATISTICS (Gelman-Rubin)
	if(rank==0) Gelman_Rubin(Nchain);

	/****************************************************************************/
	// RETURNING
	free_char_matrix(paramnames,Nparams,200);
	free_double_vector(param_read);
	free_double_vector(par_step_read);

	ierr = MPI_Finalize();

	return error_mode;
}
