#ifndef LINALG
#define LINALG

extern void multiply_squared_matrices(double **M1, double **M2, double **RES, int N);
extern void LineariseGSLMatrix(gsl_matrix * m, double * vec, int Nrows, int Ncols);

#endif
