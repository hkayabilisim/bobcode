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
          "[--kmax number_of_vawes]\n"
          "[--perturb]\n");
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
  double mass[70000];
  char inifile[100],accfile[100],potfile[100] ;
  double perturbx = 0.0;
  double perturby = 0.0;
  double perturbz = 0.0;
  bool perturb = false;
  
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
    } else if (strcmp(argv[i],"--kmax") == 0) {
      kmax = atoi(argv[i+1]);
    } else if (strcmp(argv[i],"--perturb") == 0) {
      perturb = true;
    }
  }
  ff->kLimUserSpecified = kmax;

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
    sscanf(line, "%lf%lf%lf%lf%lf", &q[i], r[i], r[i]+1, r[i]+2,&mass[i]);
  }
  fclose(ifile);
  
  if (perturb) {
    srand(time(0));
    rand(); // skip first random number which seems not random on Mac OS
    perturbx = (2.0 * rand()/(double)RAND_MAX - 1.0 )* edge[0][0];
    perturby = (2.0 * rand()/(double)RAND_MAX - 1.0 )* edge[1][1];
    perturbz = (2.0 * rand()/(double)RAND_MAX - 1.0 )* edge[2][2];
    for (int i = 0 ; i < N ; i++) {
      r[i][0] += perturbx ;
      r[i][1] += perturby ;
      r[i][2] += perturbz ;
    }
  }
  
  o.e = (double *)malloc((L+1)*sizeof(double));
 
  msm4g_tic();
  FF_build(ff, N, edge);
  double time_build = msm4g_toc();
  
  msm4g_tic();
  energy = FF_energy(ff, N, F, r, q, NULL);
  double time_energy = msm4g_toc();
  FF_get_topGridDim(ff,M);
  
  double time_manual_sum = 0;
  printf("%-30s : %10.8f\n","time_partcl2partcl",ff->time_partcl2partcl);
  time_manual_sum += ff->time_partcl2partcl;
  //double time_grid2grid_total = 0;
  for (int l = 1 ; l <= ff->maxLevel ; l++) {
    //time_grid2grid_total += ff->time_grid2grid[l];
    printf("time_grid2grid_atlevel%d        : %10.8f\n",l,ff->time_grid2grid[l]);
    printf("time_padding_atlevel%d          : %10.8f\n",l,ff->time_padding[l]);
    time_manual_sum += ff->time_grid2grid[l];
    time_manual_sum += ff->time_padding[l];
  }
  //printf("%-30s : %10.8f\n","time_grid2grid_total",time_grid2grid_total);
  //double time_restriction_total = 0;
  for (int l = 2 ; l <= ff->maxLevel ; l++) {
    //time_restriction_total += ff->time_restriction[l];
    printf("time_restrictionfrom%dto%d       : %10.8f\n",l-1,l,
           ff->time_restriction[l]);
    time_manual_sum += ff->time_restriction[l];
  }
  //printf("%-30s : %10.8f\n","time_restriction_total",time_restriction_total);
  
  //double time_prolongation_total = 0;
  for (int l = ff->maxLevel-1 ; l >= 1 ; l--) {
    //time_prolongation_total += ff->time_prolongation[l];
    printf("time_prolongationfrom%dto%d      : %10.8f\n",l+1,l,
           ff->time_prolongation[l]);
    time_manual_sum += ff->time_prolongation[l];
  }
  //printf("%-30s : %10.8f\n","time_prolongation_total",time_prolongation_total);
  
  printf("%-30s : %10.8f\n","time_anterpolation",ff->time_anterpolation);
  printf("%-30s : %10.8f\n","time_interpolation",ff->time_interpolation);
  time_manual_sum += ff->time_anterpolation + ff->time_interpolation;
  time_manual_sum += time_build;

  printf("%-30s : %10.8f\n","time_build",time_build);
  //printf("%-30s : %10.8f\n","time_energy",time_energy);
  printf("%-30s : %10.8f\n","time_other",time_build+time_energy-time_manual_sum);
  printf("%-30s : %10.8f\n","time_total",time_build+time_energy);
  printf("%-30s : %s\n", "data",argv[1]);
  printf("%-30s : %d\n", "NumberOfLevels",FF_get_maxLevel(ff));
  printf("%-30s : %d\n", "NumberOfParticles",N);
  printf("%-30s : %-10.5f\n","Perturbationx",perturbx);
  printf("%-30s : %-10.5f\n","Perturbationy",perturby);
  printf("%-30s : %-10.5f\n","Perturbationz",perturbz);
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
    double rms_top = 0.0;
    double rms_bottom = 0.0;
    for (int i = 0; i < N; i++){
      double acci2 = acc[i][0]*acc[i][0] + acc[i][1]*acc[i][1] + acc[i][2]*acc[i][2];
      double acci = sqrt(acci2);
      max_acc = fmax(max_acc, acci);
      double errx = acc[i][0] + F[i][0]/q[i],
      erry = acc[i][1] + F[i][1]/q[i],
      errz = acc[i][2] + F[i][2]/q[i];
      double err2 = errx*errx + erry*erry + errz*errz;
      double err = sqrt(err2);
      max_acc_err = fmax(max_acc_err, err);
      rms_top += err2/mass[i];
      rms_bottom += acci2/mass[i];
    }
    double rms = sqrt(rms_top/rms_bottom);
    printf("%-30s : %25.16e\n", "forceerror",max_acc_err/max_acc);
    printf("%-30s : %25.16e\n", "forcermserror",rms);
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
