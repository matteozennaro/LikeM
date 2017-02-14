#ifndef MCMC
#define MCMC

extern double randrange(double minval, double maxval);
extern void updated_covariance(gsl_matrix *C, int Nrow, int Ncols, double **par ,int current_step,int Nparams);
extern double box_muller(double mu, double sigma);
extern void compute_proposed(double *theta_1,double*theta_0,double **par,gsl_matrix *current_covariance,gsl_matrix *step_cov,gsl_matrix *mixing, gsl_vector *diag_step_cov,int Nparams,int i_step);
extern void Gelman_Rubin(int Nchain, int Nparams, char chaindir[], char root_name[], char **paramnames);

#endif
