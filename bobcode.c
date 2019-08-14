#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include "forcefield.h"

struct O o = {"csCl"}; // struct for optional output

void usage() {
  fprintf(stderr,"Usage: msm dataFile "
          "[--nbar relativeCutoff] "
          "[--nu accuracyOrder] "
          "[-M grid-spacing] \n"
          "[-L numberOfLevels] "
          "[--tol-dir direct_tolerance] "
          "[--tol-rec reciprocal_tolerance] \n"
          "[--enable-ewald-splitting] "
          "[--klim number_of_vawes]\n");
  exit(1);
}
int main(int argc, char **argv){
  if (argc == 1) {
    usage();
  }
 
  
  FF *ff = FF_new();
  int M[3] = {0, 0, 0};
  double energy;
  double edge[3][3];
  double r[70000][3];
  double F[70000][3];
  double acc[70000][3];
  double q[70000];
  char inifile[100],accfile[100],potfile[100] ;
  
  double ewaldSplitting = false ;
  int  L=-1, kmax = -1 ;
  for (int i = 0 ; i < argc ;i++) {
    if (strcmp(argv[i],"--nbar") == 0) {
      double nbar = atof(argv[i+1]);
      FF_set_relCutoff(ff, nbar);
    } else if (strcmp(argv[i],"--nu") == 0) {
      int nu = atoi(argv[i+1]);
      FF_set_orderAcc(ff, nu);
    } else if (strcmp(argv[i],"-M") == 0) {
      int Mtop = atoi(argv[i+1]);
      M[0] = Mtop; M[1] = Mtop; M[2] = Mtop;
      FF_set_topGridDim(ff, M);
    } else if (strcmp(argv[i],"-L") == 0) {
      int L = atoi(argv[i+1]);
      FF_set_maxLevel(ff, L);
    } else if (strcmp(argv[i],"--tol-dir") == 0) {
      double tol_dir = atof(argv[i+1]);
      FF_set_tolDir(ff, tol_dir);
    } else if (strcmp(argv[i],"--tol-rec") == 0) {
      double tol_rec = atof(argv[i+1]);
      FF_set_tolRec(ff, tol_rec);
    } else if (strcmp(argv[i],"--enable-ewald-splitting") == 0) {
      ewaldSplitting = true;
    } else if (strcmp(argv[i],"--kmax") == 0) {
      kmax = atoi(argv[i+1]);
    }
  }
  ff->kLimUserSpecified = kmax;
  FF_set_ewaldSplitting(ff,ewaldSplitting);


  sprintf(inifile,"%s.ini",argv[1]);
  sprintf(accfile,"%s.acc",argv[1]);
  sprintf(potfile,"%s.pot",argv[1]);
  FILE *ifile = fopen(inifile, "r");
  
  // read it in by line number and character position
  char line[180];
  
  fgets(line, sizeof line, ifile);
  sscanf(line,"%lf%lf%lf%lf%lf%lf%lf%lf%lf",
     edge[0],edge[0]+1,edge[0]+2,
     edge[1],edge[1]+1,edge[1]+2,
     edge[2],edge[2]+1,edge[2]+2);
  fgets(line, sizeof line, ifile);
  int N; sscanf(line,"%d", &N);
  for (int i = 0; i < N; i++) {
    fgets(line, sizeof line, ifile);
    sscanf(line, "%lf%lf%lf%lf", &q[i], r[i], r[i]+1, r[i]+2);
  }
  fclose(ifile);
  o.e = (double *)malloc((L+1)*sizeof(double));
 
  msm4g_tic();
  FF_build(ff, N, edge);
  double time_build = msm4g_toc();
  
  msm4g_tic();
  energy = FF_energy(ff, N, F, r, q, NULL);
  double time_energy = msm4g_toc();
  FF_get_topGridDim(ff,M);
  
  printf("%-30s : %10.8f\n","time_direct",ff->time_partcl2partcl);
  printf("%-30s : %10.8f\n","time_build",time_build);
  printf("%-30s : %10.8f\n","time_energy",time_energy);
  printf("%-30s : %10.8f\n","time_total",time_build+time_energy);
  printf("%-30s : %s\n", "data",argv[1]);
  printf("%-30s : %d\n", "NumberOfLevels",FF_get_maxLevel(ff));
  printf("%-30s : %f\n", "nbar",FF_get_relCutoff(ff));
  printf("%-30s : %d\n", "nu",FF_get_orderAcc(ff));
  printf("%-30s : %f\n", "cutoff",FF_get_cutoff(ff));
  printf("%-30s : %5.2f %5.2f %5.2f\n", "Edge row1",
           edge[0][0],edge[0][1],edge[0][2]);
  printf("%-30s : %5.2f %5.2f %5.2f\n", "Edge row2",
           edge[1][0],edge[1][1],edge[1][2]);
  printf("%-30s : %5.2f %5.2f %5.2f\n", "Edge row3",
           edge[2][0],edge[2][1],edge[2][2]);
  printf("%-30s : %d\n", "TopLevelMx",M[0]);
  printf("%-30s : %d\n", "TopLevelMy",M[1]);
  printf("%-30s : %d\n", "TopLevelMz",M[2]);
  printf("%-30s : %10.3e\n", "tol_dir",FF_get_tolDir(ff));
  printf("%-30s : %10.3e\n", "tol_rec",FF_get_tolRec(ff));
  printf("%-30s : %.16f\n", "beta",ff->beta);
  printf("%-30s : %.16f\n", "kmax",ff->kmax);
  printf("%-30s : %3d\n","klimx",ff->kLim[0]);
  printf("%-30s : %3d\n","klimy",ff->kLim[1]);
  printf("%-30s : %3d\n","klimz",ff->kLim[2]);
  printf("%-30s : %3d\n","kLimUserSpecified",ff->kLimUserSpecified);
  printf("%-30s : %3d\n","effectiveklim_x",ff->kLimUserSpecified > -1 ? ff->kLimUserSpecified : ff->kLim[0]);
  printf("%-30s : %3d\n","effectiveklim_y",ff->kLimUserSpecified > -1 ? ff->kLimUserSpecified : ff->kLim[1]);
  printf("%-30s : %3d\n","effectiveklim_z",ff->kLimUserSpecified > -1 ? ff->kLimUserSpecified : ff->kLim[2]);
  printf("%-30s : %d\n", "ewaldSplitting",ff->ewaldSplitting);
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
