#ifndef TINKER
#define TINKER

struct HMFDATA {
	int N;
	double * BinMin;
	double * BinMax;
	double * BinMean;
	double * HMF;
	double * HMFErr;
	double **C;
	double **InvC;
	double Mmin;
	double Mmax;
	double hSim;
	double VolSim;
	double Omega_m;
	double rho_bg;
};

struct STab {
	int N;
	double Mmin;
  double Mmax;
  double * M;
  double * Sigma;
};

extern void read_histo(char file[],int N,double *min,double *max,double *mean,double *mf,double *err);
extern void read_hmf_cov(char covf[], int N, double **C);
extern void compute_inverse_covariance(double **C, double **InvC, int N);
extern void FillSigmaTab(char file[], struct STab * SigmaTab);
extern void FreeSigmaTab(struct STab * SigmaTab);
extern void FreeHMFdata(struct HMFDATA * D);
extern double pkinterp(int n, double *x, double *y, double x0);
extern double W2tophat(double k, double R);
extern double func(double lnk, double lnP, double R);
extern double Sigma(double M, int k_num, double *k, double *Pk, double rho_bg);
extern double SigmaInterp(double M, struct STab * SigmaTab);
extern double ComputeNu(double M, struct STab * S);
extern void ComputeSigmaTab(char pkf[], char sigmaf[], int N, double Mmin, double Mmax, double rho_bg);
extern void InitialiseTinker(char ini[], struct HMFDATA * D, struct STab * S);
extern void QuitTinker(struct HMFDATA * HMFdata,struct STab * SigmaTab);
extern double chi2_Tinker(double *params, int Nparams, struct HMFDATA * HMFdata,struct STab * SigmaTab);

#endif
