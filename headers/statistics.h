#ifndef STATISTICS
#define STATISTICS

extern double W2(double k, double R);
extern double func1(double lnk, double lnP, double r, double R);
extern double func2(double lnk, double lnP, double R);
extern double ETA(int k_num, double *k_o, double *Pk_o, double r, double R);
extern double SIGMA2(int k_num, double *k_o, double *Pk_o, double R);
extern double XI(int k_num, double *k_o, double *Pk_o, double r, double R);

#endif
