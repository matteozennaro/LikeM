#ifndef FILEMANIP
#define FILEMANIP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

extern void read_matrix(gsl_matrix *C, int Nrow, int Ncols, char file[]);
extern int count_lines(char file[]);
extern int count_header_lines(char file[]);
extern int count_number_of_columns(char file[], int number_of_header_lines);
extern void output_for_getdist(int index, char chaindir[], char root_name[], char **paramnames, char **paramlatex, int Nparams);
extern void adjust_path(char *str);
extern void trim(char *str);
extern void create_external(char chaindir[], char root_name[], char **paramnames, char **paramlatex, int Nparams);

#endif
