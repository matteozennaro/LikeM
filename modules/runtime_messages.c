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

void frame(char str[])
{
  int add_newline=0;
  int i=0;
  while (str[i]!='\0') i++; if(str[i-1]!='\n') add_newline=1;
  printf("--------------------------------------------------------------------------\n");
  printf("%s",str);
  if (add_newline==1) printf("\n");
  printf("--------------------------------------------------------------------------\n");
  fflush(stdout);
}

int check_status(int mystatus)
{
  MPI_Barrier(MPI_COMM_WORLD);

  int everything_is_good = 1;

  if(rank!=0) MPI_Send(&mystatus,1,MPI_INT,0,1,MPI_COMM_WORLD);
  else
  {
    if(mystatus!=SUCCESS) everything_is_good = 0;
    int i;
    for (i=1;i<size;i++)
    {
      MPI_Recv(&mystatus,1,MPI_INT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      if(mystatus!=SUCCESS) everything_is_good = 0;
    }
  }

  int i;
  if(rank==0) for (i=1;i<size;i++) MPI_Send(&everything_is_good,1,MPI_INT,i,2,MPI_COMM_WORLD);
  else for (i=1;i<size;i++) if(rank==i) MPI_Recv(&everything_is_good,1,MPI_INT,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  return everything_is_good;
}

void mpi_exit(int mystatus)
{
  if(rank==0)
  {
    char str[LONG_STR];

    switch (mystatus) {
      case IO_FAILURE: {sprintf(str,"I/O error"); break;}
      case MEM_FAILURE: {sprintf(str,"Memory handling error"); break;}
      case NUM_FAILURE: {sprintf(str,"Numerical error"); break;}
      case FAILURE: {sprintf(str,"Generic failure"); break;}
      default: sprintf(str,"Unknown error occurred!");
    }

    frame(str);
  }

  MPI_Abort(MPI_COMM_WORLD,mystatus);

  MPI_Barrier(MPI_COMM_WORLD);
}

void read_err(char paramname[])
{
  char error[LONG_STR];
  sprintf(error,"Error! Parameter %s hasn't been specified\n",paramname);
  frame(error);
  status = IO_FAILURE;
}

void print_help()
{
  if(rank==0)
  {
    printf("Usage of LikeM:\n"
          "./LikeM params.ini\n\n"
          "params.ini must contain:\n\n"
          "Nparams = \n"
          "Nstep = \n"
          "root_name = \n"
          "chain_dir = \n"
          "append =   #defalt is F\n"
          "#for each parameter a line like this:\n"
          "# par[0...N] = name guess initial_jump min max latex_name\n"
          "par0  =  alpha  0.15  0.01  0.0  120.0  $ \\alpha $ \n"
          "par1  =  beta   2.31  0.02 -3.0    3.0  $ \\beta (x) $ \n"
          "# etc ...\n\n"
          "See the provided params.ini for examples.\n");
  }
}
