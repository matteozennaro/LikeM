#ifndef COLE05PK
#define COLE05PK

struct DataPk {
  int Nk;
  double *k;
  double *pk;
  double *err;
};

extern void InitialiseCole05DataPk(char parfile[], struct DataPk * data, struct DataPk * lin);
extern void QuitCole05Pk(struct DataPk * data, struct DataPk * lin);
extern double chi2_Cole05Pk(double *params, int Nparams, struct DataPk * data, struct DataPk * lin);

#endif
