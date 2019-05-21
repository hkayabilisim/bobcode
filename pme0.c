// file pme0.c
// compile with -lfftw3

// invoke methods in following order:
//   FF_new, FF_set_<parm>, FF_build
// then invoke following in any order:
//   FF_get_<parm>, FF_energy, FF_rebuild
// finally invoke FF_delete
#include <assert.h>  //:::
#include <stdio.h>  //: might remove later
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "pme.h"

typedef struct Vector {double x, y, z;} Vector;
typedef struct Matrix {double xx, xy, xz, yx, yy, yz, zx, zy, zz;} Matrix;
typedef struct Triple {int x, y, z;} Triple;

FF *FF_new(void){
  FF *ff = (FF *)calloc(1, sizeof(FF));
  return ff;}

// for each computational parameter
// have a get method FF_set_... that returns the parameter
// Each set method is optional
// No method parameter should be changed after a call to build
void FF_set_cutoff(FF *ff, double cutoff){
  ff->cutoff = cutoff;}
// ::: YET TO DO: PARAMETER CHECKING :::
// ::: e.g., the tolDir and tolRec
void FF_set_orderAcc(FF *ff, int orderAcc ){
  // orderAcc: even integer >= 4
  ff->orderAcc = orderAcc;}
void FF_set_topGridDim(FF *ff, int topGridDim[3]){
  *(Triple *)ff->topGridDim = *(Triple *)topGridDim;}
void FF_set_tolDir(FF *ff, double tolDir){
  ff->tolDir = tolDir;}
void FF_set_tolRec(FF *ff, double tolRec){
  ff->tolRec = tolRec;}

// helper functions:
static double invert(Matrix *A);
static void omegap(FF *ff);
void FF_build(FF *ff, int N, double edges[3][3]){
  // N: number of particles
  // edges: columns are vectors defining parallelpiped
  ff->N = N;
  Matrix Ai = *(Matrix *)edges;
  double detA = invert(&Ai);

  // set default values for unspecified method parameters
  if (! ff->cutoff) ff->cutoff = 8.;
  if (! ff->orderAcc) ff->orderAcc = 6;
  if (! ff->tolDir){
    if (ff->tolRec) ff->tolDir = ff->tolRec;
    else ff->tolDir = 1.e-5;}
  if (! ff->tolRec)
    ff->tolRec = ff->tolDir;
  
	// calculate beta
  // erfc(beta a_0)/a_0 = ff->tolDir/h_0  % EPBD05 omit 1/h_star
  double a_0 = ff->cutoff;
  double pi = 4.*atan(1.);
  double hstar = pow(detA/(double)N, 1./3.);
	double const_ = ff->tolDir*a_0/hstar;
  double beta = 0.; // beta = beta*a_0 until after iteration
  double res = erfc(beta) - const_;
  double oldres = 2.*res;
  while (fabs(res) < fabs(oldres)){  // Newton-Raphson
    double dres = -2./sqrt(pi)*exp(-beta*beta);
    beta -= res/dres;
    oldres = res;
    res = erfc(beta) - const_;}
  beta = beta/a_0;
  ff->beta = beta;

	if (! ff->topGridDim[0]) {
    double gridsize[2][9] = {
      0.70, 0.77, 0.83, 0.91, 1.00, 1.11, 1.21, 1.33, 1.48,
      1.05, 1.18, 1.29, 1.43, 1.54, 1.67, 1.82, 2.00, 2.22};
    int i = (ff->orderAcc - 4)/2;
    if (i > 1){
      printf("default grid dimension unavailable for order = %d\n",
             ff->orderAcc);
      exit(1);}
    int j = (int)round((ff->cutoff - 6.)/0.5);
    if (j < 0 || j > 8){
      printf("default grid dimension unavailable for cutoff = %f\n",
             ff->cutoff);
      exit(1);}
		double h_ref = gridsize[i][j];
		// calculate beta_ref
		double tolDir_ref = pow(10., -5)*hstar;
		const_ = tolDir_ref*a_0/hstar;
		beta = 0.; // beta = beta*a_0 until after iteration
		res = erfc(beta) - const_;
		oldres = 2.*res;
		while (fabs(res) < fabs(oldres)){  // Newton-Raphson
			double dres = -2./sqrt(pi)*exp(-beta*beta);
			beta -= res/dres;
			oldres = res;
			res = erfc(beta) - const_;}
		beta = beta/a_0;
		double beta_ref = beta;
		beta = ff->beta;
		double h = pow(ff->tolDir/tolDir_ref, 1./(double)ff->orderAcc)
			*(beta_ref/beta)*h_ref;
		Matrix A = *(Matrix *)edges;
		double ax = sqrt(A.xx*A.xx + A.xy*A.xy + A.xz*A.xz);
		double ay = sqrt(A.yx*A.yx + A.yy*A.yy + A.yz*A.yz);
		double az = sqrt(A.zx*A.zx + A.zy*A.zy + A.zz*A.zz);
    ff->topGridDim[0] = (int)ceil(ax/h);
    ff->topGridDim[1] = (int)ceil(ay/h);
    ff->topGridDim[2] = (int)ceil(az/h);}

  // set kmax
  const_ = sqrt(pi)*ff->tolRec/(2.*beta*hstar);
  double kmax = 0.; // kmax = pi*kmax/beta until after iteration
  res = erfc(kmax) - const_;
  oldres = 2.*res;
  while (fabs(res) < fabs(oldres)){  // Newton-Raphson
    double dres = -2./sqrt(pi)*exp(-kmax*kmax);
    kmax -= res/dres;
    oldres = res;
    res = erfc(kmax) - const_;}
  kmax = beta*kmax/pi;
  ff->kmax = kmax;
  double asx = sqrt(Ai.xx*Ai.xx + Ai.xy*Ai.xy + Ai.xz*Ai.xz);
  double asy = sqrt(Ai.yx*Ai.yx + Ai.yy*Ai.yy + Ai.yz*Ai.yz);
  double asz = sqrt(Ai.zx*Ai.zx + Ai.zy*Ai.zy + Ai.zz*Ai.zz);
  ff->kLim[0] = (int)(kmax/asx);
  ff->kLim[1] = (int)(kmax/asy);
  ff->kLim[2] = (int)(kmax/asz);
  
  // build first ord/2 pieces of B-splines Q_ord(t)
  // piece_i(t) = q_i0 + q_i1*(t - i) + ... + q_{i,ord-1}*(t - i)^{ord-1}
  int nu = ff->orderAcc, nknots = nu/2;
  ff->Q = (double *)calloc(nknots*nu, sizeof(double));
  // for degrees i = 0, 1, ..., nu-1
  // compute values of Q_{i+1} at knots j = 0, 1, ..., nknots-1
  double *Q = (double *)calloc(nu*nknots, sizeof(double));
  // Q_1(0) = 1, Q_1(j) = 0 for j > 0
  Q[0] = 1.;
  double *Qip1 = Q;
  // Q_{i+1}(j) = (j*Q_i(j) + (i+1-j)*Q_i(j-1))/i
  for (int i = 1; i < nu; i++){
    double *Qi = Qip1;
    Qip1 = Qi + nknots;
    Qip1[0] = 0.;
    for (int j = 1; j < nknots; j++)
      Qip1[j] =
        ((double)j*Qi[j] + (double)(i + 1 - j)*Qi[j-1])/(double)i;}
  // compute k-th derivative of Q_{k+1}, ..., Q_nu divided by k!
  // Q_{i+1}^(k)(j+)/k! = (Q_i^(k-1)(j+)/(k-1)! - Q_i^(k-1)(j-1+)/(k-1)!)/k
  for (int k = 1; k < nu; k++){
    // (k-1)st derivative of Q_{i+k} is in row i
    // replace with k-th derivative of Q{i+k+1}
    for (int i = 0; i < nu - k; i++){
      double *Qi = Q + i*nknots;
      for (int j = nknots - 1; j > 0; j--)
        Qi[j] = (Qi[j] - Qi[j-1])/(double)k;
      Qi[0] = Qi[0]/(double)k;}}
  for (int j = 0; j < nknots; j++)
    for (int k = 0; k < nu; k++)
      ff->Q[j*nu + k] = Q[(nu - 1 - k)*nknots + j];
  free(Q);

  // build cL
  omegap(ff);
  
  // build grid-to-grid stencil
  // assume DFT used
  // allocate memory for khat
  // special treatment for level l = maxLevel
  int di = ff->topGridDim[0], dj = ff->topGridDim[1], dk = ff->topGridDim[2];
  ff->khat = (double *)calloc(di*dj*dk, sizeof(double));
  ff->fftw_in = (fftw_complex *)fftw_malloc(di*dj*dk*sizeof(fftw_complex));
  ff->forward = fftw_plan_dft_3d(di, dj, dk, ff->fftw_in, ff->fftw_in,
                                FFTW_FORWARD, FFTW_ESTIMATE);
  ff->backward = fftw_plan_dft_3d(di, dj, dk, ff->fftw_in, ff->fftw_in,
                                FFTW_BACKWARD, FFTW_ESTIMATE);

  // compute stencil ff->khat
  FF_rebuild(ff, edges);
}

// for each computational parameter
// also have a get method FF_set_... that returns the parameter
double FF_get_cutoff(FF*ff) {
  return ff->cutoff;}
int FF_get_orderAcc(FF*ff) {
  return ff->orderAcc;}
void FF_get_topGridDim(FF*ff, int topGridDim[3]) {
  *(Triple *)topGridDim = *(Triple *)ff->topGridDim;}
double FF_get_tolDir(FF*ff) {
  return ff->tolDir;}
double FF_get_tolRec(FF*ff) {
  return ff->tolRec;}

double FF_get_errEst(FF *ff, int N, double *charge){
  // calculate C_{nu-1}
  int nu = ff->orderAcc;
  if (nu > 10){
    printf("No error estimate available for order %d.\n", nu);
    return nan("");}
  int index = (nu - 4)/2;
  double M[4] = {9., 825., 130095., 34096545.};
  double cbarp[4] = {1./6., 1./30., 1./140., 1./630.};
  double k[4] = {0.39189561, 0.150829428, 0.049632967, 0.013520855};
  double C = 4./3.*M[index]*cbarp[index]*k[index];
  // calculate hmin
  Matrix Ai = *(Matrix *)ff->A;
  double detA = invert(&Ai);
  double asx = sqrt(Ai.xx*Ai.xx + Ai.xy*Ai.xy + Ai.xz*Ai.xz);
  double asy = sqrt(Ai.yx*Ai.yx + Ai.yy*Ai.yy + Ai.yz*Ai.yz);
  double asz = sqrt(Ai.zx*Ai.zx + Ai.zy*Ai.zy + Ai.zz*Ai.zz);
  Triple M1 = *(Triple *)ff->topGridDim;
  Vector h
    = {1./(asx*(double)M1.x), 1./(asy*(double)M1.y), 1./(asz*(double)M1.z)};
  double hmin = fmin(fmin(h.x, h.y), h.z);
  // complete the calculation
  double Q2 = 0;
  for (int i = 0; i < N; i++) Q2 += charge[i]*charge[i];
  double a0 = ff->cutoff;
  return Q2*C*pow(hmin, nu-2)/(pow(detA, 1./3.)*sqrt((double)N)*pow(a0, nu-1));
}

// The rebuild method is called every time the periodic cell changes.
// It initializes edges and calculate the grid2grid stencils.
// helper functions:
static void dALp1(FF *ff, Triple gd, double kh[], double detA);
void FF_rebuild(FF *ff, double edges[3][3]) {
  *(Matrix *)ff->A = *(Matrix *)edges;
  Matrix A = *(Matrix *)ff->A;
  // calculate coeffs for const part of energy
  int nu = ff->orderAcc;
  double a_0 = ff->cutoff;
  double beta = ff->beta;
  double pi = 4.*atan(1.);
  Matrix Ai = *(Matrix *)edges;
  double detA = invert(&Ai);
  *(Matrix *)ff->Ai = Ai;
  Vector as = {sqrt(Ai.xx*Ai.xx + Ai.yx*Ai.yx + Ai.zx*Ai.zx),
               sqrt(Ai.xy*Ai.xy + Ai.yy*Ai.yy + Ai.zy*Ai.zy),
               sqrt(Ai.xz*Ai.xz + Ai.yz*Ai.yz + Ai.zz*Ai.zz)};
  double xlim = ceil(a_0*as.x - 1.);
  double ylim = ceil(a_0*as.y - 1.);
  double zlim = ceil(a_0*as.z - 1.);
  ff->coeff1 = 0.;  // coeff of 1/2(sum_i q_i^2)
  for (double x = -xlim; x <= xlim; x++)
    for (double y = -ylim; y <= ylim; y++)
      for (double z = -zlim; z <= zlim; z++){ // could do just z >= 0
        double rx = A.xx*x + A.xy*y + A.xz*z;
        double ry = A.yx*x + A.yy*y + A.yz*z;
        double rz = A.zx*x + A.zy*y + A.zz*z;
        double r2 = rx*rx + ry*ry + rz*rz;
        if (r2 >= a_0*a_0) continue;
        if (r2 != 0.){
          double r1 = sqrt(r2);
          ff->coeff1 += erfc(beta*r1)/r1;}
        else ff->coeff1 -= 2.*beta/sqrt(pi);
        //***
        //***
        }
  // compute f->coeff2
  ff->coeff2 = pi/(beta*beta*detA);  // coeff of 1/2(sum_i q_i)^2

  // build grid-to-grid stencil for levels L, L + 1
  // ff->khat[L] = d^{L+1} + DFT of khat^L
  Triple gd = *(Triple *)ff->topGridDim;
  double *kh = ff->khat;
  // add in d^2(A)
  // d^2(A)_n = sum_k chi(k) c'(k) exp(2 pi i k . H_L n)
  dALp1(ff, gd, kh, detA);
}

void FF_delete(FF *ff) {
  free(ff->Q);
  free(ff->khat);
  for (int alpha = 0; alpha < 3; alpha++) free(ff->cL[alpha]);
  fftw_destroy_plan(ff->forward);
  fftw_destroy_plan(ff->backward);
  fftw_free(ff->fftw_in);
  free(ff);
}

//helper functions

static double invert(Matrix *A) {  // invert A
  // returns |det A|
  double (*a)[3] = (double (*)[3])A;
  // Matrix Ai; double (*ai)[3] = (double (*)[3])&Ai;
  double ai[3][3];
  // compute adjugate
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      ai[i][j] = a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]
                - a[(j+2)%3][(i+1)%3]*a[(j+1)%3][(i+2)%3];
  // A adj(A) = (det A) I
  double detA = a[0][0]*ai[0][0] + a[0][1]*ai[1][0] + a[0][2]*ai[2][0];
  for (int j = 0; j < 9; j++) a[0][j] = ai[0][j]/detA;
  return fabs(detA);}

static void omegap(FF *ff){
  // construct ff->cL
  int nu = ff->orderAcc;
  // nu positive even integer
  int nuby2 = nu/2;

  // Compute Phi(n), n = 0, 1,..., nuby2-1
  double *Phi = (double *)calloc(nuby2, sizeof(double));
  for(int i = 1; i < nuby2; i++)  Phi[nuby2-i] = ff->Q[i*nu];
  for (int k = 0; k < nu; k++) Phi[0] += ff->Q[(nuby2-1)*nu + k];

  int M[3] = {ff->topGridDim[0], ff->topGridDim[1], ff->topGridDim[2]};
  double pi = 4.*atan(1.);
  for (int alpha = 0; alpha < 3; alpha++){
    int Ma = M[alpha];
    double *c = (double *)malloc((Ma/2 + 1)*sizeof(double));
    for (int k = 0; k <= Ma/2; k++){
      c[k] = Phi[0];
      double theta = 2.*pi*(double)k/(double)Ma;
      for (int n = 1; n <= nu/2 - 1; n++)
        c[k] += 2.*Phi[n]*cos((double)n*theta);
      c[k] *= c[k];
      c[k] = 1./c[k];
    }
    ff->cL[alpha] = c;
  }
 free(Phi);
}

static void dALp1(FF *ff, Triple gd, double kh[], double detA){
  // add d^{L+1}(A) to kh
  Matrix Ai = *(Matrix *)ff-> Ai;
  Triple kLim = *(Triple *)ff->kLim;
  double pi = 4.*atan(1.);
  double pidetA = pi*fabs(detA);
  double pi2beta2 = pow(pi/ff->beta, 2);
  // loop on vec k
  // for kx = 0, 1, -1, ..., kLim.x, -kLim.x
  for (int kx = - kLim.x ; kx <= kLim.x; kx++){
    int kx1 = (kx % gd.x + gd.x) % gd.x;
    int kx0 = (kx1 + gd.x/2) % gd.x - gd.x/2;
    double cLx = ff->cL[0][abs(kx0)];
    for (int ky = - kLim.y; ky <= kLim.y; ky++){
      int ky1 = (ky % gd.y + gd.y) % gd.y;
      int ky0 = (ky1 + gd.y/2) % gd.y - gd.y/2;
      double cLxy = cLx*ff->cL[1][abs(ky0)];
      for (int kz = - kLim.z; kz <= kLim.z; kz++){
        int kz1 = (kz % gd.z + gd.z) % gd.z;
        int kz0 = (kz1 + gd.z/2) % gd.z - gd.z/2;
        double cLxyz = cLxy*ff->cL[2][abs(kz0)];
        double fkx = (double)kx, fky = (double)ky, fkz = (double)kz;
        double Aikx = Ai.xx*fkx + Ai.yx*fky + Ai.zx*fkz,
          Aiky = Ai.xy*fkx + Ai.yy*fky + Ai.zy*fkz,
          Aikz = Ai.xz*fkx + Ai.yz*fky + Ai.zz*fkz;
        double k2 = Aikx*Aikx + Aiky*Aiky + Aikz*Aikz;
        double chi_cL = k2 == 0 ? 0 : exp(- pi2beta2*k2)/(pidetA*k2)*cLxyz;
        kh[(kx1*gd.y + ky1)*gd.z + kz1]
          += chi_cL*gd.x*gd.y*gd.z;
      }
    }
  }
}
