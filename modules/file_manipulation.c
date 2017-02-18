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
#include <runtime_messages.h>
#include <vectors.h>
#include <read_ini_file.h>

void read_matrix(gsl_matrix *C, int Nrow, int Ncols, char file[])
{
  FILE *f = fopen(file,"r");
  if (f==NULL)
  {
    frame("You need to specify a covariance matrix!");
    status = IO_FAILURE;
  }
  else
  {
    int i,j;
    double val;
    for (i=0; i<Nrow; i++)
    {
      for (j=0; j<Ncols; j++)
      {
        fscanf(f,"%lf",&val);
        gsl_matrix_set(C,i,j,val);
      }
    }
    fclose(f);
  }
}

int count_lines(char file[])
{
  int l=0;
  FILE *f = fopen(file,"r");
  if (f==NULL)
  {
    char error[LONG_STR];
    sprintf(error,"Error opening file %s\n",file);
    frame(error);
    status = IO_FAILURE;
  }
  else
  {
    char test;
    while (fscanf(f,"%c",&test)!=EOF) if(test=='\n') l++;
    fclose(f);
  }
  return l;
}

int count_header_lines(char file[])
{
  int hdl=0;
  char buf[LONG_STR];
  char * token;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL) exit(-1);
      token = strtok(buf," \t");
      if (token[0]=='#' || token[0]=='/' || token[0]=='!') hdl++;
      else break;
    }
    fclose(f);
  }
  else
  {
    char error[LONG_STR];
    sprintf(error,"Problem loading file %s!\n",file);
    frame(error);
    status=IO_FAILURE;
  }
  return hdl;
}

int count_number_of_columns(char file[], int number_of_header_lines)
{
  int ncol=0;
  int line=0;
  char buf[LONG_STR];
  char * token;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    while(line<number_of_header_lines)
    {
      if(fgets(buf,sizeof(buf),f)==NULL) exit(-1);
      line++;
    }
    if(fgets(buf,sizeof(buf),f)==NULL) exit(-1);
    token = strtok(buf," \t\n");
    if(token!=NULL)
    {
      while (token!=NULL)
      {
        ncol++;
        token = strtok(NULL," \t\n");
      }
    }
    fclose(f);
  }
  else
  {
    char error[LONG_STR];
    sprintf(error,"Problem loading file %s!\n",file);
    frame(error);
    status = IO_FAILURE;
  }
  return ncol;
}


void output_for_getdist(int index, char chaindir[], char root_name[], char **paramnames, char **paramlatex, int Nparams)
{
  char outfile[STD_STR];
  sprintf(outfile,"%s%s_%i_raw.txt",chaindir,root_name,index);
  FILE *f = fopen(outfile,"r");
  if (f==NULL)
  {
    char error[LONG_STR];
    sprintf(error,"Error opening file %s",outfile);
    frame(error);
    status = IO_FAILURE;
  }
  else
  {
    int nl = count_lines(outfile);
    if(status!=SUCCESS) return;
    int n_burn_in = (int)(0.05*nl);

    double **par = allocate_double_matrix(Nparams,nl);
    if(status!=SUCCESS) return;

    double par_now[Nparams];

    double chi2_now;
    double *chi2_vec = allocate_double_vector(nl);
    if(status!=SUCCESS) return;

    int *N = allocate_int_vector(nl);
    if(status!=SUCCESS) return;

    int Nidependent=0;

    int i,j,indexpar,indexline;
    double line_diff;
    int found_matching_line = 0;

    for (i=0; i<nl; i++) N[i] = 0;

    int k=0;

    for (i=0; i<nl; i++)
    {
      fscanf(f,"%lf",&chi2_now);
      for (j=0; j<Nparams; j++) fscanf(f,"%lf",&par_now[j]);

      if (i>n_burn_in)
      {
        for (indexline=0; indexline<Nidependent; indexline++)
        {
          line_diff = 0.0;
          for (indexpar=0; indexpar<Nparams; indexpar++)
          line_diff += (par[indexpar][indexline]-par_now[indexpar]);
          if (fabs(line_diff) < 1.0e-9) {N[indexline]++; found_matching_line++;}
        }

        //printf("%i\n",found_matching_line);

        if (found_matching_line == 0.0)
        {
          for (indexpar=0; indexpar<Nparams; indexpar++)
          {
            par[indexpar][k] = par_now[indexpar];
          }
          chi2_vec[k] = chi2_now;
          N[k]++;
          Nidependent++;
          k++;
        }
        else {found_matching_line=0;}
      }
    }
    fclose(f);

    double chi2min=chi2_vec[0];
    for (i=0; i<Nidependent; i++)
    {
      if (chi2_vec[i] < chi2min) chi2min = chi2_vec[i];
    }

    double ParMin[Nparams];
    double ParMax[Nparams];
    for (indexpar=0; indexpar<Nparams; indexpar++)
    {
      ParMin[indexpar] = par[indexpar][0];
      ParMax[indexpar] = par[indexpar][0];
    }

    sprintf(outfile,"%s%s_%i.txt",chaindir,root_name,index);
    FILE * out = fopen(outfile,"w");
    if(out==NULL)
    {
      printf("Error creating file %s\n",outfile);
      status = IO_FAILURE;
      return;
    }
    for (i=0; i<Nidependent; i++)
    {
      fprintf(out,"%i\t%.10e\t",N[i],0.5*(chi2_vec[i]-chi2min));
      for (indexpar=0; indexpar<Nparams; indexpar++)
      {
        fprintf(out,"%.10e\t",par[indexpar][i]);

        if (par[indexpar][i]>ParMax[indexpar])
          ParMax[indexpar]=par[indexpar][i];
        if (par[indexpar][i]<ParMin[indexpar])
          ParMin[indexpar]=par[indexpar][i];
      }
      fprintf(out,"\n");
    }
    fclose(out);

    sprintf(outfile,"%s%s.paramnames",chaindir,root_name);
    out = fopen(outfile,"w");
    for(indexpar=0;indexpar<Nparams;indexpar++)
    {
      fprintf(out,"%s\t%s\n",paramnames[indexpar],paramlatex[indexpar]);
    }
    fclose(out);

    sprintf(outfile,"%s%s.ranges",chaindir,root_name);
    out = fopen(outfile,"w");
    for (indexpar=0;indexpar<Nparams;indexpar++)
    {
      fprintf(out,"%s  %e  %e\n",
      paramnames[indexpar],
      ParMin[indexpar]-0.1*ParMin[indexpar],
      ParMax[indexpar]+0.1*ParMax[indexpar]);
    }
    fclose(out);

    // char command[2000];
    // sprintf(command,"mv %s_%i.txt chains",root_name,index);
    // system(command);
    // sprintf(command,"rm %s_%i_raw.txt* ",root_name,index);
    // system(command);

    deallocate_double_matrix(par,Nparams,nl);
    deallocate_int_vector(N);
    deallocate_double_vector(chi2_vec);
  }
}

void adjust_path(char *str)
{
  int i=0;
  while(str[i]!='\0') i++;
  if(str[i-1]!='/')
  {
    str[i] = '/';
    str[i+1] = '\0';
  }
}

void trim(char *str)
{
  int i=0;
  char *tmp;
  while(str[i]!='\0') i++;
  while(str[i-1]==' ' || str[i-1]=='\t') {i--; str[i]='\0';}
  i=0;
  while(str[i]==' ' || str[i]=='\t') i++;
  tmp = str+i;
  sprintf(str,"%s",tmp);
}

void create_external(char chaindir[], char root_name[], char **paramnames,
  char **paramlatex, int Nparams)
{
  char getdistfile[STD_STR];
  char textablefile[STD_STR];
  int i;

  sprintf(getdistfile,"external/getdist_%s.ini",root_name);
  sprintf(textablefile,"external/tex_table_%s.ini",root_name);

  FILE *f;
  f = fopen(getdistfile,"w");
  if(f==NULL)
  {
    printf("Error creating file %s\n",getdistfile);
    status = IO_FAILURE;
    return;
  }
  else
  {
    fprintf(f,
      "no_tests = F\n"
      "file_root = %s/%s\n"
      "out_root = dist_%s\n"
      "out_dir = %s\n"
      "plot_data_dir = ./plot_data/\n"
      "chain_num = -1\n"
      "first_chain =\n"
      "exclude_chain =\n"
      "ignore_rows = 0.3\n"
      "#include defaults settings for kernel densitiy estimates etc, can also be specified in this file if you want to override\n"
      "DEFAULT(/home/matteo/CosmoCodes/CosmoMC/python/getdist/analysis_defaults.ini)\n"
      "samples_are_chains = T\n"
      "no_plots = T\n"
      "plot_2D_param = 0\n"
      "plot_2D_num = 0\n"
      "plot1 =\n"
      "plot2 =\n"
      "num_3D_plots = 0\n"
      "3D_plot1 =\n"
      "marker[nrun] = 0\n"
      "thin_factor = 0\n"
      "make_single_samples = F\n"
      "single_thin = 4\n"
      "limits[r02]= 0 N\n"
      "limits[r10]= 0 N\n"
      "all_limits =\n"
      "force_twotail = F\n"
      "PCA_num = 0\n"
      "PCA_normparam =\n"
      "PCA_params =\n"
      "PCA_func   = LLL\n"
      "cool = 1\n",
      chaindir,root_name,
      root_name,
      chain);
    fclose(f);
  }

  f = fopen(textablefile,"w");
  if(f==NULL)
  {
    printf("Error creating file %s\n",textablefile);
    status = IO_FAILURE;
    return;
  }
  else
  {
    fprintf(f,
      "Ncosmo = 1\n"
      "cosmo_dir = %s\n"
      "cosmo_names = dist_%s\n"
      "cosmo_tex = %s\n"
      "Nparams = %i\n",
      chaindir,
      root_name,
      root_name,
      Nparams);
    fprintf(f,
      "param_names = ");
    for(i=0;i<Nparams;i++) fprintf(f,"%s ",paramnames[i]); fprintf(f,"\n");
    fprintf(f,
      "param_tex = ");
    for(i=0;i<Nparams;i++) fprintf(f,"%s ",paramlatex[i]); fprintf(f,"\n");
    fprintf(f,
      "out_file = table_%s.tex\n"
      "# Optional parameters:\n"
      "overwrite = T #default is F\n"
      "create_tex_doc = T  #default is F\n"
      "split_table =  T #default is F\n"
      "print_bestfits = T #default is F\n",
      root_name);
}
