#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>
#include <ctype.h>

#include <extern_global.h>
#include <runtime_messages.h>
#include <file_manipulation.h>

int read_double_from_file(char file[], char paramname[], double *val)
{
  int reading_success=0;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[LONG_STR];
    char *token;

    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL)
      {
        return reading_success;
      }
      token = strtok(buf," \t=\n");
      if (token!=NULL)
      {
        if(token[0]!='#')
        {
          while (token!=NULL)
          {
            if (strcmp(token,paramname)==0)
            {
              token = strtok(NULL," \t=\n");
              *val = atof(token);
              reading_success = 1;
              if (verb>0) printf("  %s = %lf\n",paramname,*val);
              break;
            }
            else
            {
              token = strtok(NULL," \t=\n");
            }
          }
        }
      }
      if (reading_success==1) break;
    }
    fclose(f);
  }
  else
  {
    char error[LONG_STR];
    sprintf(error,"Error! Unable to open param file %s\n"
          "Make sure you have specified a valid param file\n",file);
    frame(error);
    status = IO_FAILURE;
  }
  return reading_success;
}

int read_int_from_file(char file[], char paramname[], int *val)
{
  int reading_success=0;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[LONG_STR];
    char *token;

    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL)
      {
        return reading_success;
      }
      token = strtok(buf," \t=\n");
      if (token!=NULL)
      {
        if(token[0]!='#')
        {
          while (token!=NULL)
          {
            if (strcmp(token,paramname)==0)
            {
              token = strtok(NULL," \t=\n");
              *val = atoi(token);
              reading_success = 1;
              if (verb>0) printf("  %s = %i\n",paramname,*val);
              break;
            }
            else
            {
              token = strtok(NULL," \t=\n");
            }
          }
        }
      }
      if (reading_success==1) break;
    }
    fclose(f);
  }
  else
  {
    char error[LONG_STR];
    sprintf(error,"Error! Unable to open param file %s\n"
          "Make sure you have specified a valid param file\n",file);
    frame(error);
    status = IO_FAILURE;
  }
  return reading_success;
}

int read_bool_from_file(char file[], char paramname[], char *truth)
{
  int reading_success=0;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[LONG_STR];
    char *token;

    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL)
      {
        return reading_success;
      }
      token = strtok(buf," \t=\n");
      if (token!=NULL)
      {
        if(token[0]!='#')
        {
          while (token!=NULL)
          {
            if (strcmp(token,paramname)==0)
            {
              token = strtok(NULL," \t=\n");
              if (token[0]=='T' || token[0]=='t')
                *truth = 'T';
              else
                *truth = 'F';
              reading_success = 1;
              if (verb>0) printf("  %s = %c\n",paramname,*truth);
              break;
            }
            else
            {
              token = strtok(NULL," \t=\n");
            }
          }
        }
      }
      if (reading_success==1) break;
    }
    fclose(f);
  }
  else
  {
    char error[LONG_STR];
    sprintf(error,"Error! Unable to open param file %s\n"
          "Make sure you have specified a valid param file\n",file);
    frame(error);
    status = IO_FAILURE;
  }
  return reading_success;
}

int read_string_from_file(char file[], char paramname[], char *str)
{
  int reading_success=0;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[LONG_STR];
    char *token;

    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL)
      {
        return reading_success;
      }
      token = strtok(buf," \t=\n");
      if (token!=NULL)
      {
        if(token[0]!='#')
        {
          while (token!=NULL)
          {
            if (strcmp(token,paramname)==0)
            {
              token = strtok(NULL," \t=\n");
              sprintf(str,"%s",token);
              reading_success = 1;
              if (verb>0) printf("  %s = %s\n",paramname,str);
              break;
            }
            else
            {
              token = strtok(NULL," \t=\n");
            }
          }
        }
      }
      if (reading_success==1) break;
    }
    fclose(f);
  }
  else
  {
    char error[LONG_STR];
    sprintf(error,"Error! Unable to open param file %s\n"
          "Make sure you have specified a valid param file\n",file);
    frame(error);
    status = IO_FAILURE;
  }
  return reading_success;
}

int read_param_from_file(char file[],char paramname[], char * name, char *latex,
  double *p, double * jump, double *parmin, double *parmax)
{
  int reading_success=0;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[LONG_STR];
    char *token;

    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL)
      {
        return reading_success;
      }
      token = strtok(buf," \t=\n");
      if (token!=NULL)
      {
        if(token[0]!='#')
        {
          while (token!=NULL)
          {
            if (strcmp(token,paramname)==0)
            {
              token = strtok(NULL," \t=\n");
              if(token!=NULL) sprintf(name,"%s",token);
              else {status = IO_FAILURE; return reading_success;}

              token = strtok(NULL," \t=\n");
              if(token!=NULL) *p = atof(token);
              else {status = IO_FAILURE; return reading_success;}

              token = strtok(NULL," \t=\n");
              if(token!=NULL) *jump = atof(token);
              else {status = IO_FAILURE; return reading_success;}

              token = strtok(NULL," \t=\n");
              if(token!=NULL) *parmin = atof(token);
              else {status = IO_FAILURE; return reading_success;}

              token = strtok(NULL," \t=\n");
              if(token!=NULL) *parmax = atof(token);
              else {status = IO_FAILURE; return reading_success;}

              token = strtok(NULL,"\n");
              if(token!=NULL) {sprintf(latex,"%s",token); trim(latex);}
              else {status = IO_FAILURE; return reading_success;}

              reading_success = 1;
              if (verb>0) printf("  %s (%s) = %lf %lf\n",paramname,name,*p,*jump);
              break;
            }
            else
            {
              token = strtok(NULL," \t=\n");
            }
          }
        }
      }
      if (reading_success==1) break;
    }
    fclose(f);
  }
  else
  {
    char error[LONG_STR];
    sprintf(error,"Error! Unable to open param file %s\n"
          "Make sure you have specified a valid param file\n",file);
    frame(error);
    status = IO_FAILURE;
  }
  return reading_success;
}
