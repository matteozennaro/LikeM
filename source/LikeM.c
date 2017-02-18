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
#include <file_manipulation.h>
#include <linalg.h>
#include <mcmc_functions.h>
#include <read_ini_file.h>
#include <runtime_messages.h>
#include <statistics.h>
#include <vectors.h>

#include <TinkerHMF.h>
#include <Cole05bias.h>

int main (int argc, char * argv[])
{
	status = SUCCESS;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank==0)
	{
		// STARTING EVERYTHING...
		if (argc == 2)
		{
			if (strcmp(argv[1],"--help")==0)
			{
				print_help();
				MPI_Finalize();
				return 0;
			}
			else
			{
				printf("\n*************************************************************\n");
				printf("\n                Welcome in LikeM!\n\n");
				printf("*************************************************************\n");

				getcwd(workdir,sizeof(workdir));
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
			frame("Error! You must specify an input parameter file\n"
						"For further info on the usage of mcmc, please type\n"
						"./mcmc --help\n");
			mpi_exit(FAILURE);
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

  if(!read_int_from_file(argv[1],"verb",&verb)) verb = 1;
  if (rank!=0) verb = 0;

  char append;
  int Nparams;
  int Nchain = size;
  int Nstep;
  if(!read_bool_from_file(argv[1],"append",&append)) append = 'F';
  if(!read_int_from_file(argv[1],"Nparams",&Nparams)) read_err("Nparams");
  if(!read_int_from_file(argv[1],"Nstep",&Nstep)) read_err("Nstep");
	if(!check_status(status)) mpi_exit(status);

  char *root_name = allocate_char_vector(STD_STR);
  char *chaindir = allocate_char_vector(LONG_STR);
	if(!check_status(status)) mpi_exit(status);

  if(!read_string_from_file(argv[1],"root_name",root_name)) read_err("root_name");
	if(!check_status(status)) mpi_exit(status);
	if(!read_string_from_file(argv[1],"chain_dir",chaindir)) sprintf(chaindir,"");
  else adjust_path(chaindir);

  double *param_read = allocate_double_vector(Nparams);
	double *par_step_read = allocate_double_vector(Nparams);
	double *par_min_read = allocate_double_vector(Nparams);
  double *par_max_read = allocate_double_vector(Nparams);
	char **paramnames = allocate_char_matrix(Nparams,STD_STR);
  char **paramlatex = allocate_char_matrix(Nparams,STD_STR);
	if(!check_status(status)) mpi_exit(status);

  int i;
  char par_tag[STD_STR];
  for (i=0; i<Nparams; i++)
  {
    sprintf(par_tag,"par%i",i);
    if(!read_param_from_file(argv[1],par_tag,paramnames[i],paramlatex[i],&param_read[i],&par_step_read[i],&par_min_read[i],&par_max_read[i]))
    {
      char error[LONG_STR];
      sprintf(error,"Error! Reading of parameter par%i failed\n",i);
			frame(error);
      mpi_exit(IO_FAILURE);
    }
  }

	// now let's understand which likelihoods are required
	if(!read_bool_from_file(argv[1],"use_tinker",&use_tinker)) use_tinker='F';
	else if(use_tinker=='T') LikeNum++;

	if(!read_bool_from_file(argv[1],"use_cole05",&use_cole05)) use_cole05='F';
	else if(use_cole05=='T') LikeNum++;

	if(LikeNum == 0)
	{
		char error[LONG_STR];
		sprintf(error,"Error! In your param file you didn't specify any\n"
						"likelihood to be used. This is a fatal error.");
		if(rank==0) frame(error);
		status = FAILURE;
	}

	if(!check_status(status)) mpi_exit(status);
	MPI_Barrier(MPI_COMM_WORLD);

	NU = 2.38/sqrt(Nparams);

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

  gsl_matrix *current_covariance;

  double **par;

  char outfile[STD_STR];
  time_t start,stop;

  for(i_par=0;i_par<Nparams;i_par++)
  {
    param[i_par] = param_read[i_par];
    par_step[i_par] = par_step_read[i_par];
  }

  /****************************************************************************/
  // CHECKS
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
		if (param[i_par] < par_min_read[i_par] || param[i_par] > par_max_read[i_par])
		goto redoparam;
	}

	#ifdef DEBUG
		Nstep = 3;
		for (i=0; i<Nparams; i++)
		{
			printf("%.10e	%.10e\n",param[i],par_step[i]);
		}
	#endif

  /****************************************************************************/
  // INITIALIZE RELEVANT LIKELIHOODS
	struct HMFDATA HMFdata;
	struct STab SigmaTab;
	if (use_tinker=='T')
	{
		InitialiseTinker(argv[1],&HMFdata,&SigmaTab);
		if(!check_status(status)) mpi_exit(status);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	struct BiasData data;
	if (use_cole05=='T')
	{
		InitialiseCole05Data(argv[1],&data);
		if(!check_status(status)) mpi_exit(status);
		MPI_Barrier(MPI_COMM_WORLD);
	}

  /****************************************************************************/
	// METROPOLIS - HASTINGS
	start = time(NULL);

	N_accepted = N_tot = 0;
	sprintf(outfile,"%s%s_%i_raw.txt",chaindir,root_name,rank+1);

	i_step=0;
	printf("Starting chain %i!\n",rank+1);
	fflush(stdout);

	for (i_par=0; i_par<Nparams; i_par++) theta_0[i_par] = param[i_par];

	chi2_0 = 0.0;
	if(use_tinker=='T') chi2_0 += chi2_Tinker(theta_0,Nparams,&HMFdata,&SigmaTab);
	if(!check_status(status)) mpi_exit(status);

	if(use_cole05=='T') chi2_0 += chi2_Cole05(theta_0,Nparams,&data);
	if(!check_status(status)) mpi_exit(status);

	if (append == 'T') // NEED TO UPDATE THIS
	{
		// retrieve last point from previously computed chain
		int nl = count_lines(outfile);
		if(!check_status(status)) mpi_exit(status);

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

		printf("Last line in chain %i\n",rank+1);
		for (i_par=0; i_par<Nparams; i_par++) printf("%e  ",theta_0[i_par]);
		printf("\n");
	}

	current_covariance = gsl_matrix_alloc(Nparams,Nparams);
	par = allocate_double_matrix(Nstep,Nparams);
	if(!check_status(status)) mpi_exit(status);

	FILE * out;
	if (append=='T') out=fopen(outfile,"a");
	else out=fopen(outfile,"w");
  if(out==NULL)
  {
    char error[LONG_STR];
    sprintf(error,"Error! Unable to open file %s\n",outfile);
    frame(error);
		status = IO_FAILURE;
  }
	if(!check_status(status)) mpi_exit(status);

	for (i_par=0; i_par<Nparams; i_par++) par[0][i_par] = theta_0[i_par];
	for (i_par=0; i_par<Nparams; i_par++) gsl_matrix_set(current_covariance,i_par,i_par,par_step[i_par]*par_step[i_par]);

	int ii,jj;

	gsl_matrix * mixing = gsl_matrix_alloc(Nparams,Nparams);
	gsl_matrix * step_cov = gsl_matrix_alloc(Nparams,Nparams);
	gsl_vector * diag_step_cov = gsl_vector_alloc(Nparams);

	for (i_step=0; i_step<Nstep; i_step++)
	{
		respect_priors:
		if (i_step>999)
		{
			compute_proposed(theta_1,theta_0,par,current_covariance,step_cov,mixing,diag_step_cov,Nparams,i_step);
		}
		else
		{
			for (i_par=0; i_par<Nparams; i_par++) theta_1[i_par] = theta_0[i_par] + par_step[i_par]*randrange(-1.0,1.0);
		}
		for (i_par=0; i_par<Nparams; i_par++)
		if(theta_1[i_par]<par_min_read[i_par] || theta_1[i_par]>par_max_read[i_par])
		goto respect_priors;


		chi2_1 = 0.0;
		if(use_tinker=='T') chi2_1 += chi2_Tinker(theta_1,Nparams,&HMFdata,&SigmaTab);
		if(!check_status(status)) mpi_exit(status);

		if(use_cole05=='T') chi2_1 += chi2_Cole05(theta_1,Nparams,&data);
		if(!check_status(status)) mpi_exit(status);

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

  // Wait for all the chains to finish
  MPI_Barrier(MPI_COMM_WORLD);

	int proc = size;
  for(i=1;i<proc;i++)
  {
    if(i==rank)
    {
      int PassInfo[3];
      PassInfo[0] = N_accepted;
      PassInfo[1] = N_tot;
      PassInfo[2] = (int) (stop-start);
      MPI_Send(PassInfo,3,MPI_INT,0,0,MPI_COMM_WORLD);

      double *LinearisedCovariance = allocate_double_vector(Nparams*Nparams);
      LineariseGSLMatrix(current_covariance,LinearisedCovariance,Nparams,Nparams);
      MPI_Send(LinearisedCovariance,Nparams*Nparams,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
      deallocate_double_vector(LinearisedCovariance);
    }
  }

  if(rank==0)
  {
    printf("\n");
  	printf("\nCovariance from chain %i:\n\n",rank+1);
  	for (ii=0;ii<Nparams;ii++)
  	{
  		for(jj=0;jj<Nparams;jj++) printf("%.3e  ",gsl_matrix_get(current_covariance,ii,jj));
  		printf("\n");
  	}
  	printf("\n");

  	printf("Ch %i, Number of accepted steps = %i\n"
  				 "Ch %i, Acceptance rate = %lf\n"
  				 "Ch %i, Total elapsed time = %ld s.\n",
  				 rank+1,N_accepted,
  				 rank+1,(double)N_accepted/N_tot,
  				 rank+1,stop-start);

    int RcvInfo[3];
    double *RcvCov = allocate_double_vector(Nparams*Nparams);
    for(i=1;i<proc;i++)
    {
      MPI_Recv(RcvInfo,3,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(RcvCov,Nparams*Nparams,MPI_DOUBLE,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

      printf("\n");
  	   printf("\nCovariance from chain %i:\n\n",i+1);
    	for (ii=0;ii<Nparams;ii++)
    	{
    		for(jj=0;jj<Nparams;jj++) printf("%.3e  ",RcvCov[ii*Nparams+jj]);
    		printf("\n");
    	}
    	printf("\n");

    	printf("Ch %i, Number of accepted steps = %i\n"
    				 "Ch %i, Acceptance rate = %lf\n"
    				 "Ch %i, Total elapsed time = %i s.\n",
    				 rank+1,RcvInfo[0],
    				 rank+1,(double)RcvInfo[0]/RcvInfo[1],
    				 rank+1,RcvInfo[2]);
    }
    deallocate_double_vector(RcvCov);
  }


	/****************************************************************************/
	// REARRANGE OUTPUT FILE FOR USING getdist
	output_for_getdist(rank+1,chaindir,root_name,paramnames,paramlatex,Nparams);
	if(!check_status(status)) mpi_exit(status);

	if(rank==0) create_external(chaindir,root_name,paramnames,paramlatex,Nparams);
	if(!check_status(status)) mpi_exit(status);

	/****************************************************************************/
	// PREPARING TO QUIT
	gsl_matrix_free(current_covariance);
	gsl_matrix_free(mixing);
	gsl_matrix_free(step_cov);
	gsl_vector_free(diag_step_cov);
	deallocate_double_matrix(par,Nstep,Nparams);

  // QUIT RELEVANT LIKELIHOODS
	if (use_tinker=='T') QuitTinker(&HMFdata,&SigmaTab);
	if (use_cole05=='T') QuitCole05(&data);

	// COMPUTE CONVERGENCE STATISTICS (Gelman-Rubin)
	if(rank==0) Gelman_Rubin(Nchain,Nparams,chaindir,root_name,paramnames);
	if(!check_status(status)) mpi_exit(status);

	/****************************************************************************/
	// RETURNING
	deallocate_char_matrix(paramnames,Nparams,STD_STR);
	deallocate_double_vector(param_read);
	deallocate_double_vector(par_step_read);

	MPI_Finalize();

	return 0;
}
