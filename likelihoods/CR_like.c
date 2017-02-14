#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

extern double ETA(int k_num, double *k, double *Pk, double r, double R);

extern double *ytab;
extern double *FFtab;
extern double *GGtab;
extern double ytab_max;
extern double ytab_min;
extern double F_inf,F_0;
extern double G_inf;
extern double ytabstep;
extern int ntab;

extern int error_mode;
extern char classdir[500];
extern char dir[500];

int count_lines(char file[])
{
  int l=0;
  char test='a';
  FILE *f = fopen(file,"r");
  if (f == NULL)
  {
    printf("Error opening file %s\n",file);
    error_mode = -1;
  }
  else
  {
    while(fscanf(f,"%c",&test)!=EOF) if (test=='\n') l++;
    fclose(f);
  }
  return l;
}

void read_GG_FF_tabs()
{
  char dir[500];
  sprintf(dir,"tabulated_functions/FF_GG_func_tab.dat");
  FILE *gft = fopen(dir,"r");
  if (gft==NULL)
  {
    printf("\n\tOoops! A folder named \'tabulated_functions\' was expected to be found\n"
    "\there %s\n"
    "\tand it was expected to contain 1 file named \'FF_GG_func_tab.dat\' but \n"
    "\tsomething went wrong.\n\n",dir);
    error_mode=-1;
    return ;
  }

  int i; int l=0;
  char test='a';
  while(fscanf(gft,"%c",&test)!=EOF) if (test=='\n') l++;
  fseek(gft,0,SEEK_SET);

  ytab = (double*)malloc(l*sizeof(double));
  GGtab = (double*)malloc(l*sizeof(double));
  FFtab = (double*)malloc(l*sizeof(double));
  if (ytab==NULL || GGtab == NULL || FFtab == NULL)
  {
    printf("Bad memory allocation!\n");
    error_mode = -1;
    return ;
  }

  for (i=0; i<l; i++)
  {
    fscanf(gft,"%lf %lf %lf",&ytab[i],&FFtab[i],&GGtab[i]);
  }
  fclose(gft);
  ytab_min = ytab[0];
  ytab_max = ytab[l-1];
  ntab = l;
  ytabstep = log(ytab[1])-log(ytab[0]);

  double x=0.;
  double xmin=0.;
  double xmax=20.;
  double xstep=(xmax-xmin)/(100000.-1.);

  F_inf=0.;
  for(x=xmin; x<=xmax; x+=xstep)
  {
    F_inf += (x*x/(1.+exp(x)))*xstep;
  }

  G_inf=F_inf;

  F_0=0.;
  for(x=xmin; x<=xmax; x+=xstep)
  {
    F_0 += (x*x*x/(1.+exp(x)))*xstep;
  }
}

double FF(double Y)
{
  int i=0;
  double m,q;
  if (Y >= ytab_max)
  {
    return Y*F_inf;
  }
  else if (Y < ytab_min)
  {
    return F_0;
  }
  else
  {
    i = 0;
    i = floor((log(Y)-log(ytab_min))/(ytabstep));
    m = (FFtab[i]-FFtab[i+1])/(ytab[i]-ytab[i+1]);
    q = FFtab[i] - m*ytab[i];
    return (m*Y+q);
  }
}

double GG(double Y)
{
  int i=0;
  double m,q;
  if (Y >= ytab_max)
  {
    return G_inf;
  }
  else if (Y < ytab_min)
  {
    i=0;
    m=(GGtab[i]-0.)/(ytab[i]-0.);
    q = 0.;
    return (m*Y+q);
  }
  else
  {
    i=floor((log(Y)-log(ytab_min))/(ytabstep));
    m = (GGtab[i]-GGtab[i+1])/(ytab[i]-ytab[i+1]);
    q = GGtab[i] - m*ytab[i];
    return (m*Y+q);
  }
}

double ONE2(double z, double mnu, double h)
{
  if (mnu==0.0) return 0.0;
  else
  {
    double Gamma_nu = 0.71611;
    double A = 1.0/(1.0+z);
    double N_nu = 3.0;
    double Tcmb_0 = 2.7255;
    double Kb = 8.617342e-5;
    double Gnu4 = Gamma_nu*Gamma_nu*Gamma_nu*Gamma_nu;
    double pi4 = M_PI*M_PI*M_PI*M_PI;
    double a4 = A*A*A*A;
    double y = (mnu*A)/(Gamma_nu*Kb*N_nu*Tcmb_0);
    return ((15.0*Gnu4*N_nu*(2.469e-05/(h*h)))/(a4*pi4))*FF(y);
  }
}

double func_1_3w(double z, double mnu)
{
  if (mnu==0.0) return 1.0;
  else
  {
    double Gamma_nu = 0.71611;
    double A = 1.0/(1.0+z);
    double N_nu = 3.0;
    double Tcmb_0 = 2.7255;
    double Kb = 8.617342e-5;
    double y = (mnu*A)/(Gamma_nu*Kb*N_nu*Tcmb_0);
    return y*(GG(y)/FF(y));
  }
}

void create_class_ini(double *theta, int nz, double z[nz], int chain_ind)
{
  // MAKE THINGS HUMAN READABLE!!!
  double H,OB0h2,OC0h2,mnu,tau,As,ns;
  H = 100.0*theta[0];
  OB0h2 = theta[1];
  OC0h2 = theta[2];
  mnu = theta[3];
  tau = theta[4];
  As = theta[5];
  ns = theta[6];

  double kmax = 50;

  double N_ur,N_ncdm;
  if (mnu == 0.0)
  {
    N_ur = 3.046;
    N_ncdm = 0.0;
  }
  else
  {
    N_ur = 0.00641;
    N_ncdm = 3.0;
  }

  char filename[200];
  sprintf(filename,"power_%i.ini",chain_ind);

  FILE *out = fopen(filename,"w");

  fprintf(out,
  "H0 = %.8lf                      \n"
  "T_cmb = 2.7255                      \n"
  "omega_b = %.8lf                     \n"
  "N_ur = %.8lf                      \n"
  "omega_cdm = %.8lf                     \n"
  "Omega_dcdmdr = 0.0                    \n"
  "Gamma_dcdm = 0.0                    \n"
  "N_ncdm = %lf                      \n"
  "m_ncdm = %4lf,%4lf,%4lf                   \n"
  "Omega_ncdm =                        \n"
  "T_ncdm = 0.71611,0.71611,0.71611                \n"
  "deg_ncdm = 1,1,1                      \n"
  "Omega_k = 0.                      \n"
  "#Omega_Lambda = 0.7                     \n"
  "                        \n"
  "attractor_ic_scf = yes                    \n"
  "#scf_parameters = [scf_lambda, scf_alpha, scf_A, scf_B, phi, phi_prime]       \n"
  "scf_parameters = 10.0, 0.0, 0.0, 0.0, 100.0, 0.0            \n"
  "scf_tuning_index = 0                    \n"
  "                        \n"
  "YHe = BBN                       \n"
  "recombination = RECFAST                   \n"
  "reio_parametrization = reio_camb                \n"
  "tau_reio = %.8lf                    \n"
  "reionization_exponent = 1.5                   \n"
  "reionization_width = 0.5                  \n"
  "helium_fullreio_redshift = 3.5                  \n"
  "helium_fullreio_width = 0.5                   \n"
  "binned_reio_num = 3                     \n"
  "binned_reio_z = 8,12,16                   \n"
  "binned_reio_xe = 0.8,0.2,0.1                  \n"
  "binned_reio_step_sharpness = 0.3                \n"
  "annihilation = 0.                     \n"
  "annihilation_variation = 0.                   \n"
  "annihilation_z = 1000                     \n"
  "annihilation_zmax = 2500                  \n"
  "annihilation_zmin = 30                    \n"
  "annihilation_f_halo= 20                   \n"
  "annihilation_z_halo= 8                    \n"
  "on the spot = yes                     \n"
  "decay = 0.                      \n"
  "                        \n"
  "output = mPk,mTk                    \n"
  "                        \n"
  "non linear =                        \n"
  "                        \n"
  "modes = s                       \n"
  "                        \n"
  "lensing = no                      \n"
  "                        \n"
  "ic = ad                       \n"
  "gauge = synchronous                     \n"
  "                        \n"
  "P_k_ini type = analytic_Pk                  \n"
  "k_pivot = 0.05                      \n"
  "A_s = %.8e                          \n"
  "n_s = %.8lf                       \n"
  "alpha_s = 0.                      \n"
  "                        \n"
  "f_bi = 1.                       \n"
  "n_bi = 1.5                      \n"
  "f_cdi=1.                      \n"
  "f_nid=1.                      \n"
  "n_nid=2.                      \n"
  "alpha_nid= 0.01                     \n"
  "                        \n"
  "c_ad_bi = 0.5                       \n"
  "c_ad_cdi = -1.                      \n"
  "c_bi_nid = 1.                       \n"
  "                        \n"
  "r = 1.                        \n"
  "n_t = scc                       \n"
  "alpha_t = scc                       \n"
  "                        \n"
  "potential = polynomial                    \n"
  "V_0=1.e-13                      \n"
  "V_1=-1.e-14                       \n"
  "V_2=7.e-14                      \n"
  "V_3=                        \n"
  "V_4=                        \n"
  "                        \n"
  "H_0=1.e-13                      \n"
  "H_1=-1.e-14                       \n"
  "H_2=7.e-14                      \n"
  "H_3=                        \n"
  "H_4=                        \n"
  "                        \n"
  "phi_end =                       \n"
  "Vparam0 =                       \n"
  "Vparam1 =                       \n"
  "Vparam2 =                       \n"
  "Vparam3 =                       \n"
  "Vparam4 =                       \n"
  "ln_aH_ratio = 50                    \n"
  "                        \n"
  "k1=0.002                      \n"
  "k2=0.1                        \n"
  "                        \n"
  "P_{RR}^1 = 2.3e-9                     \n"
  "P_{RR}^2 = 2.3e-9                     \n"
  "P_{II}^1 = 1.e-11                     \n"
  "P_{II}^2 = 1.e-11                     \n"
  "P_{RI}^1 = -1.e-13                    \n"
  "|P_{RI}^2| = 1.e-13                     \n"
  "                        \n"
  "special_iso =                       \n"
  "                        \n"
  "command = cat external_Pk/Pk_example.dat              \n"
  "                        \n"
  "custom1 = 0.05     # In the example command: k_pivot            \n"
  "custom2 = 2.215e-9 # In the example command: A_s            \n"
  "custom3 = 0.9624   # In the example command: n_s            \n"
  "custom4 = 2e-10    # In the example (with tensors) command: A_t         \n"
  "custom5 = -0.1     # In the example (with tensors) command: n_t         \n"
  "#custom6 = 0                      \n"
  "#custom7 = 0                      \n"
  "#custom8 = 0                      \n"
  "#custom9 = 0                      \n"
  "#custom10 = 0                       \n"
  "                        \n"
  "P_k_max_h/Mpc = %.1lf                     \n"
  "z_pk = ",
  H,
  OB0h2,
  N_ur,
  OC0h2,
  N_ncdm,
  mnu/3.0,mnu/3.0,mnu/3.0,
  tau,
  As,
  ns,
  kmax);

  int i;
  for (i=0; i<nz-1; i++) fprintf(out,"%lf, ",z[i]);
  fprintf(out,"%lf\n",z[nz-1]);

  fprintf(out,
  "                        \n"
  "selection=gaussian                    \n"
  "selection_mean = 0.98,0.99,1.0,1.1,1.2                \n"
  "selection_width = 0.1                     \n"
  "non_diagonal=4                      \n"
  "                        \n"
  "dNdz_selection =                    \n"
  "dNdz_evolution =                    \n"
  "bias = 1.                       \n"
  "                        \n");

  fprintf(out,"root = %s/power_%i_\n",dir,chain_ind);

  fprintf(out,
  "headers = no                      \n"
  "format = class                      \n"
  "                        \n"
  "write background = no                     \n"
  "write thermodynamics = no                   \n"
  "write primordial = no                     \n"
  "write parameters = no                           \n"
  "write warnings =                    \n"
  "                        \n"
  "input_verbose = 0                     \n"
  "background_verbose = 0                    \n"
  "thermodynamics_verbose = 0                  \n"
  "perturbations_verbose = 0                   \n"
  "transfer_verbose = 0                    \n"
  "primordial_verbose = 0                    \n"
  "spectra_verbose = 0                     \n"
  "nonlinear_verbose = 0                     \n"
  "lensing_verbose = 0                     \n"
  "output_verbose = 0                    \n");
  fclose(out);
}

void select_PCB_only(int index, int knum, double *theta, int chain_ind)
{
  // MAKE IT HUMAN READABLE!!!!!
  double h,OB0h2,OC0h2,mnu,tau,As,ns;
  h = theta[0];
  OB0h2 = theta[1];
  OC0h2 = theta[2];
  mnu = theta[3];
  tau = theta[4];
  As = theta[5];
  ns = theta[6];

  double kpivot = 0.05;

  char file[500];
  sprintf(file,"power_%i_z%i_tk.dat",chain_ind,index);
  FILE * f = fopen(file,"r");
  if (f == NULL)
  {
    error_mode = -1;
    printf("\n\tError opening file %s\n\n",file);
    return;
  }
  else
  {
    char dummy[1000];
    fgets(dummy,sizeof(dummy),f);

    double val;
    double k,pk,tc,tb,tcb;

    sprintf(file,"power_%i_z%i_PCB.dat",chain_ind,index);
    FILE * out_pcb = fopen(file,"w");

    int i;
    if (mnu==0.0)
    {
      for (i=0; i<knum; i++)
      {
        fscanf(f,"%lf %lf %lf %lf %lf %lf",&k,&val,&tb,&tc,&val,&val);
        tcb = (OB0h2/(OB0h2+OC0h2))*tb + (OC0h2/(OB0h2+OC0h2))*tc;
        pk = As*pow((k*h)/kpivot,ns-1.0)*((2.0*M_PI*M_PI)/(k*k*k))*tcb*tcb;
        fprintf(out_pcb,"%.10e\t%.10e\n",k,pk);
      }
    }
    else
    {
      for (i=0; i<knum; i++)
      {
        fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&k,&val,&tb,&tc,&val,&val,&val,&val,&val);
        tcb = (OB0h2/(OB0h2+OC0h2))*tb + (OC0h2/(OB0h2+OC0h2))*tc;
        pk = As*pow((k*h)/kpivot,ns-1.0)*((2.0*M_PI*M_PI)/(k*k*k))*tcb*tcb;
        fprintf(out_pcb,"%.10e\t%.10e\n",k,pk);
      }
    }
    fclose(out_pcb);
    fclose(f);
  }
}

void create_cosmo_file(double * theta, double z, int chain_ind)
{
  // MAKE IT HUMAN READABLE!!!!!
  double h,OB0h2,OC0h2,mnu,tau,As,ns;
  h = theta[0];
  OB0h2 = theta[1];
  OC0h2 = theta[2];
  mnu = theta[3];
  tau = theta[4];
  As = theta[5];
  ns = theta[6];

  char file[500];
  sprintf(file,"power_%i.cosmo",chain_ind);

  FILE * f = fopen(file,"w");

  fprintf(f,"Oc0          %lf\n"
  "Ob0          %lf\n"
  "mnu          %lf\n"
  "h            %lf\n"
  "z            %lf\n",
  OC0h2/(h*h),
  OB0h2/(h*h),
  mnu,
  h,
  z);

  fclose(f);
}

void read_pk(char file[],double *k, double *P, int knum)
{
  FILE *f = fopen(file,"r");
  if (f == NULL)
  {
    printf("\n\tError opening file %s\n",file);
    error_mode = -1;
    return;
  }
  else
  {
    int i;
    for (i=0; i<knum; i++)
    {
      fscanf(f,"%lf %lf",&k[i],&P[i]);
    }
    fclose(f);
  }
}

double angular_diameter_distance(double z, double *theta)
{
  //MAKE IT HUMAN READABLE!!!!!
  double h,OB0,OC0,mnu,tau,As,ns;
  h = theta[0];
  OB0 = theta[1]/(h*h);
  OC0 = theta[2]/(h*h);
  mnu = theta[3];
  tau = theta[4];
  As = theta[5];
  ns = theta[6];

  double c = 299792.458; // [km/s]

  double Neff;
  if (mnu==0.0) Neff = 3.046;
  else Neff = 0.00641;

  double OR0 = (2.469e-05/(h*h))*(1.0+Neff*(7.0/8.0)/pow(4.0/11.0,4.0/3.0));
  double OM0 = OB0+OC0+mnu/(93.14*h*h);
  double OL0 = 1.0 - OM0 - OR0;

  double x,zz;
  int nstep = 5e5;
  double step = z/(nstep-1);
  double left,right,E;
  double D = 0.0;

  int i;
  for (i=0; i<nstep-1; i++)
  {
    x = 0.0 + i*step;
    zz = 1.0+x;
    E = sqrt(OR0*zz*zz*zz*zz + OM0*zz*zz*zz + OL0);
    left = 1/E;
    zz = 1.0+x+step;
    E = sqrt(OR0*zz*zz*zz*zz + OM0*zz*zz*zz + OL0);
    right = 1/E;
    D += 0.5 * step * (left+right);
  }
  D *= c/100.0;
  return D;
}

double ap_aplha(double *theta, double * param_fiducial, double z)
{
  // MAKE IT HUMAN READABLE!!!!!
  double h,OB0,OC0,mnu,tau,As,ns;
  h = theta[0];
  OB0 = theta[1]/(h*h);
  OC0 = theta[2]/(h*h);
  mnu = theta[3];
  tau = theta[4];
  As = theta[5];
  ns = theta[6];

  double hF,OB0F,OC0F,mnuF,tauF,AsF,nsF;
  hF = param_fiducial[0];
  OB0F = param_fiducial[1]/(hF*hF);
  OC0F = param_fiducial[2]/(hF*hF);
  mnuF = param_fiducial[3];
  tauF = param_fiducial[4];
  AsF = param_fiducial[5];
  nsF = param_fiducial[6];

  double E,EF;
  double DA,DAF;

  double zz3 = (1.0+z)*(1.0+z)*(1.0+z);
  double zz4 = zz3*(1.0+z);

  double Neff;
  if (mnu==0.0) Neff = 3.046;
  else Neff = 0.00641;

  double NeffF;
  if (mnuF==0.0) NeffF = 3.046;
  else NeffF = 0.00641;

  double OR0 = (2.469e-05/(h*h))*(1.0+Neff*(7.0/8.0)/pow(4.0/11.0,4.0/3.0));
  double OM0 = OB0+OC0+mnu/(93.14*h*h);
  double OL0 = 1.0 - OM0 - OR0;

  double OR0F = (2.469e-05/(hF*hF))*(1.0+NeffF*(7.0/8.0)/pow(4.0/11.0,4.0/3.0));
  double OM0F = OB0F+OC0F+mnuF/(93.14*hF*hF);
  double OL0F = 1.0 - OM0F - OR0F;

  E = sqrt(OR0*zz4 + (OB0+OC0)*zz3 + ONE2(z,mnu,h) + OL0);
  EF = sqrt(OR0F*zz4 + (OB0F+OC0F)*zz3 + ONE2(z,mnuF,h) + OL0F);

  DA = angular_diameter_distance(z,theta);
  DAF = angular_diameter_distance(z,param_fiducial);

  return pow((EF*DA*DA)/(E*DAF*DAF),1.0/3.0);
}



double chi2_eta(double *theta, int chain_ind)
{
  #ifdef ONLYPLANCK
    return 0.0;
    printf("ciao!\n");
  #endif

  double z[4] = {0.48551,0.74658,1.05352,1.45825};

  // double eta_fid[4] = {0.1002365169,0.1005153808,0.1017705212,0.1028776456};
  double eta_fid[4] = {0.1024613033,0.1002707455,0.1018887668,0.1015149839};
  double err_eta[4] = {0.0020042834,0.0020132220,0.0019931043,0.0019996851};
  double R_fid = 22.0;
  double n_fid = 2.1;
  double param_fiducial[7] = {
    0.67,
    0.022445,
    0.117982042,
    0.3,
    0.0925,
    2.1265e-9,
    0.96 };

  double eta_theo;
  double alpha;
  double R,r;

  double chi2 = 0.0;

  create_class_ini(theta,4,z,chain_ind);

  char command[1000];
  char filename[500];
  sprintf(command,"%s/class power_%i.ini",classdir,chain_ind);
  system(command);

  int zindex;
  sprintf(filename,"power_%i_z1_pk.dat",chain_ind);
  int knum = count_lines(filename);
  if (error_mode == -1) return 0.0;

  double k[knum];
  double Pcb[knum];

  for (zindex=0; zindex<4; zindex++)
  {
    select_PCB_only(zindex+1,knum,theta,chain_ind); if (error_mode == -1) return 0.0;
    create_cosmo_file(theta,z[zindex],chain_ind);
    sprintf(command,"halofit/apply_halofit power_%i_z%i_PCB.dat power_%i.cosmo",chain_ind,zindex+1,chain_ind);
    system(command);
    sprintf(command,"power_%i_z%i_PCB_NL.dat",chain_ind,zindex+1);
    read_pk(command,k,Pcb,knum); if (error_mode == -1) return 0.0;
    alpha = ap_aplha(theta,param_fiducial,z[zindex]);
    R = alpha*R_fid;
    r = n_fid*R;
    eta_theo = ETA(knum,k,Pcb,r,R);
    chi2 += ((eta_theo - eta_fid[zindex])/err_eta[zindex])
          * ((eta_theo - eta_fid[zindex])/err_eta[zindex]);
#ifdef DEBUG 
    printf("z = %lf\n",z[zindex]);
    printf("alpha = %lf\n",alpha);
    printf("R = %lf\n",R);
    printf("eta theo = %e\n",eta_theo);
    printf("eta meas = %e\n",eta_fid[zindex]);
#endif
  }

  for (zindex=0; zindex<4; zindex++)
  {
    sprintf(command,"rm power_%i_z%i_pk.dat",chain_ind,zindex+1);
    system(command);
    sprintf(command,"rm power_%i_z%i_tk.dat",chain_ind,zindex+1);
    system(command);
    sprintf(command,"rm power_%i_z%i_PCB.dat",chain_ind,zindex+1);
    system(command);
    sprintf(command,"rm power_%i_z%i_PCB_NL.dat",chain_ind,zindex+1);
    system(command);
  }

  #ifdef DEBUG
    printf("chi2 eta: %.10e\n",chi2);
  #endif

  return chi2;
}
