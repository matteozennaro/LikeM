#ifndef COLE05
#define COLE05

struct BiasData {
  int Nk;
  double *k;
  double *b;
  double *bsigma;
};

extern void InitialiseCole05Data(char parfile[], struct BiasData * data);
extern void QuitCole05(struct BiasData * data);
extern double chi2_Cole05(double *params, int Nparams, struct BiasData * data);

#endif
