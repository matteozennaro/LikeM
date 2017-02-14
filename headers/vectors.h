#ifndef VECTORS
#define VECTORS

extern int * allocate_int_vector(int nelems);
extern double * allocate_double_vector(int nelems);
extern char * allocate_char_vector(int nelems);
extern double ** allocate_int_matrix(int nrows, int ncols);
extern double ** allocate_double_matrix(int nrows, int ncols);
extern char ** allocate_char_matrix(int nrows, int ncols);
extern void deallocate_int_vector(int *vec);
extern void deallocate_double_vector(double *vec);
extern void deallocate_char_vector(double *vec);
extern void deallocate_int_matrix(int **m, int nrows, int ncols);
extern void deallocate_double_matrix(double **m, int nrows, int ncols);
extern void deallocate_char_matrix(char **m, int nrows, int ncols);

#endif
