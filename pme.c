#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "pme.h"

struct O o = {"csCl"}; // struct for optional output

double msm4g_tictocmanager(int push) {
  double elapsed_seconds = 0.0;
  static clock_t stack_data[100] ;
  static int stack_lastindex = 0 ;
  if (push) {
    stack_data[stack_lastindex] = clock();
    stack_lastindex = stack_lastindex + 1;
  } else {
    clock_t now = clock();
    stack_lastindex = stack_lastindex - 1;
    clock_t previous = stack_data[stack_lastindex];
    elapsed_seconds = (double)(now-previous)/CLOCKS_PER_SEC;
  }
  return elapsed_seconds;
}
void msm4g_tic() {
  msm4g_tictocmanager(1);
}

double msm4g_toc() {
  return msm4g_tictocmanager(0);
}


int main(int argc, char **argv){
  if (argc != 7) {
    printf("Usage: %s data a0 nu M tol_dir tol_rec\n",argv[0]);
    exit(1);
  }
  double edge[3][3] = {20., 0., 0., 0., 20, 0., 0., 0., 20.};
  double r[30000][3];
  double F[30000][3];
  double acc[30000][3];
  double q[30000];
  char inifile[100],accfile[100],potfile[100] ;
  sprintf(inifile,"%s.ini",argv[1]);
  sprintf(accfile,"%s.acc",argv[1]);
  sprintf(potfile,"%s.pot",argv[1]);
  FILE *ifile = fopen(inifile, "r");
  
  // read it in by line number and character position
  char line[180];
  
  fgets(line, sizeof line, ifile);
  int N; sscanf(line,"%d", &N);
  for (int i = 0; i < N; i++) {
    fgets(line, sizeof line, ifile);
    sscanf(line, "%lf%lf%lf%lf", &q[i], r[i], r[i]+1, r[i]+2);
  }
  fclose(ifile);
  
  msm4g_tic();
  FF *ff = FF_new();
  int M[3] = {0, 0, 0};
  double energy;
  double a0 = atof(argv[2]);
  int nu = atoi(argv[3]);
  int Min = atoi(argv[4]);
  double tol_dir = atof(argv[5]);
  double tol_rec = atof(argv[6]);
  if (nu != 0 ) FF_set_orderAcc(ff, nu);
  if (tol_dir != 0) FF_set_tolDir(ff, tol_dir);
  if (tol_rec != 0) FF_set_tolRec(ff, tol_rec);
  if (Min != 0 ) { 
    M[0] = M[1] = M[2] = Min; 
    FF_set_topGridDim(ff, M); }
  if (a0 != 0) FF_set_cutoff(ff, a0);
  
  FF_build(ff, N, edge);
  energy = FF_energy(ff, N, F, r, q, NULL);
  
  FF_get_topGridDim(ff,M);
  printf("%-30s : %10.8f\n","time",msm4g_toc());
  printf("%-30s : %s\n", "data",argv[1]);
  printf("%-30s : %d\n", "nu",FF_get_orderAcc(ff));
  printf("%-30s : %f\n", "beta",ff->beta);
  printf("%-30s : %d\n", "Mx",M[0]);
  printf("%-30s : %d\n", "My",M[1]);
  printf("%-30s : %d\n", "Mz",M[2]);
  printf("%-30s : %f\n", "cutoff",FF_get_cutoff(ff));
  printf("%-30s : %10.3e\n", "tol_dir",FF_get_tolDir(ff));
  printf("%-30s : %10.3e\n", "tol_rec",FF_get_tolRec(ff));
  printf("%-30s : %.16e\n", "utotal",energy);
  
  FILE *afile = fopen(accfile, "r");
  if (afile != NULL) {
    fgets(line, sizeof line, afile);
    for (int i = 0 ; i < N ; i++) {
      fgets(line, sizeof line, afile); sscanf(line, "%lf", &(acc[i][0])); }
    for (int i = 0 ; i < N ; i++) {
      fgets(line, sizeof line, afile); sscanf(line, "%lf", &(acc[i][1])); }
    for (int i = 0 ; i < N ; i++) {
      fgets(line, sizeof line, afile); sscanf(line, "%lf", &(acc[i][2])); }
    
    double max_acc = 0.;
    double max_acc_err = 0.;
    for (int i = 0; i < N; i++){
      double acci
      = sqrt(acc[i][0]*acc[i][0] + acc[i][1]*acc[i][1] +	acc[i][2]*acc[i][2]);
      max_acc = fmax(max_acc, acci);
      double errx = acc[i][0] + F[i][0]/q[i],
      erry	= acc[i][1] + F[i][1]/q[i],
      errz	= acc[i][2] + F[i][2]/q[i];
      double err = sqrt(errx*errx + erry*erry + errz*errz);
      max_acc_err = fmax(max_acc_err, err);
    }
    printf("%-30s : %25.16e\n", "forceerror",max_acc_err/max_acc);
    fclose(afile);
  }
  
  
  FILE *pfile = fopen(potfile, "r");
  
  if (pfile != NULL) {
    double energy_expected=0.0;
    fgets(line, sizeof line, pfile);
    sscanf(line,"%lf", &energy_expected);
    printf("%-30s : %25.16e\n", "poterror",fabs(energy_expected-energy)/fabs(energy_expected));
    fclose(pfile);
  }
  
  FILE *fp = fopen("bob.acc","w");
  fprintf(fp,"%d\n",N);
  for (int i=0;i<N;i++) fprintf(fp,"%-25.16f\n",-F[i][0]/q[i]);
  for (int i=0;i<N;i++) fprintf(fp,"%-25.16f\n",-F[i][1]/q[i]);
  for (int i=0;i<N;i++) fprintf(fp,"%-25.16f\n",-F[i][2]/q[i]);
  fclose(fp);
  fp = fopen("bob.pot","w");
  fprintf(fp,"%25.16e\n",energy);
  fclose(fp);
  
  FF_delete(ff);
  
}
