#ifndef NUMERICS
#define NUMERICS

extern double min_double_vec(int N, double *vec);
extern double max_double_vec(int N, double *vec);
extern double lin_interp(double x0, double *x, double *y, int dimension);
extern int check_for_nan_and_inf(double val, char whatsthis[]);

#endif
