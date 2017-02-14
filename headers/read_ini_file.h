/**   READ INI FILE ***********************************************************/
#ifndef READINI
#define READINI

extern int read_double_from_file(char file[], char paramname[], double *val);
extern int read_int_from_file(char file[], char paramname[], int *val);
extern int read_bool_from_file(char file[], char paramname[], char *truth);
extern int read_string_from_file(char file[], char paramname[], char *str);
extern int read_param_from_file(char file[],char paramname[], char * name, char *latex, double *p, double * jump, double *parmin, double *parmax);

#endif
/******************************************************************************/
