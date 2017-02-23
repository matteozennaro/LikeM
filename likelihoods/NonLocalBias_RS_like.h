#ifndef NONLOC_RS
#define NONLOC_RS

#include <gsl/gsl_spline.h>

struct DataNonLocPk {
  int Nk;
  double *k;
  double *pk;
  double *err;
};

struct LinNonLocPk {
  int Nk;
  double *k;
  double *pk;
  gsl_interp_accel *acc;
  gsl_spline *splinePk;
};

extern void QuitNonLocRSPk(struct DataNonLocPk *data, struct LinNonLocPk *lin);
extern void InitialiseNonLocRSPk(char parfile[], struct DataNonLocPk *gal, struct DataNonLocPk *data, struct LinNonLocPk *lin);
extern double chi2_NonLocRSPk(double *params, int Nparams, struct DataNonLocPk * gal, struct DataNonLocPk * data, struct LinNonLocPk * lin);

#endif
