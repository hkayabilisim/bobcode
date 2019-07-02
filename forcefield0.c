// file forcefield.c
// compile with -lfftw3 or -DNO_FFT

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
#include "forcefield.h"

typedef struct Vector {double x, y, z;} Vector;
typedef struct Matrix {double xx, xy, xz, yx, yy, yz, zx, zy, zz;} Matrix;
typedef struct Triple {int x, y, z;} Triple;

FF *FF_new(void){
  FF *ff = (FF *)calloc(1, sizeof(FF));
#ifndef NO_FFT
  ff->FFT = true;
#endif
  return ff;}

// for each computational parameter
// have a get method FF_set_... that returns the parameter
// Each set method is optional
// No method parameter should be changed after a call to build
void FF_set_relCutoff(FF *ff, double relCutoff){
  ff->relCutoff = relCutoff;}
// ::: YET TO DO: PARAMETER CHECKING :::
// ::: e.g., the tolDir and tolRec
void FF_set_orderAcc(FF *ff, int orderAcc ){
  // orderAcc: even integer >= 4
  ff->orderAcc = orderAcc;}
void FF_set_maxLevel(FF *ff, int maxLevel) {
  // level-l cutoff relative to level-l grid size
  //: currently relative to minimum box dimension
  ff->maxLevel = maxLevel;}
void FF_set_topGridDim(FF *ff, int topGridDim[3]){
  *(Triple *)ff->topGridDim = *(Triple *)topGridDim;}
void FF_set_tolDir(FF *ff, double tolDir){
  ff->tolDir = tolDir;}
void FF_set_tolRec(FF *ff, double tolRec){
  ff->tolRec = tolRec;}
void FF_set_FFT(FF *ff, bool FFT){
#ifdef NO_FFT
  if (FFT){
    printf("FFT unavailable; execution aborted\n");
    exit(1);}
#endif
  ff->FFT = FFT;}

void FF_set_altSplitting(FF *ff,bool altSplitting) {
  ff->altSplitting = altSplitting;
}

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
  if (! ff->relCutoff) ff->relCutoff = 8.;
  if (! ff->orderAcc){
    if (ff->relCutoff <= 7.6) ff->orderAcc = 4;
    else if (ff->relCutoff <= 10.8) ff->orderAcc = 6;
    else ff->orderAcc = 8;}
  if (ff->relCutoff < ff->orderAcc/2){
    ff->relCutoff = ff->orderAcc/2;
    printf("relative cutoff too small; reset to %f\n",ff->relCutoff);}
  ff->nLim = ceil(ff->relCutoff - 1.);
  double asx = sqrt(Ai.xx*Ai.xx + Ai.xy*Ai.xy + Ai.xz*Ai.xz);
  double asy = sqrt(Ai.yx*Ai.yx + Ai.yy*Ai.yy + Ai.yz*Ai.yz);
  double asz = sqrt(Ai.zx*Ai.zx + Ai.zy*Ai.zy + Ai.zz*Ai.zz);
  if (! ff->topGridDim[0]) {
    double Ls = ff->maxLevel
      ? (double)ff->maxLevel
      : ceil(1. + log2((double)N)/6.);
    double Ms = pow(8., (1. - Ls))*(double)N;
    double asprod = asx * asy * asz;
    double asmin = fmin(asx, fmin(asy, asz));
    double Msmin = pow((double)(ff->nLim + 1)*asmin, 3)/asprod;
    Ms = fmax(Ms, Msmin);
    double numerator = pow(Ms*asprod, 1./3.);
    ff->topGridDim[0] = (int)ceil(numerator/asx);
    ff->topGridDim[1] = (int)ceil(numerator/asy);
    ff->topGridDim[2] = (int)ceil(numerator/asz);}
  double aL = ff->relCutoff/fmax(ff->topGridDim[0]*asx,
              fmax(ff->topGridDim[1]*asy,ff->topGridDim[2]*asz));
  printf("refl:%f asx:%f ff->topGridDim[0]:%d\n",ff->relCutoff,asx,ff->topGridDim[0]);
  if (! ff->maxLevel){
    int M = ff->topGridDim[0]*ff->topGridDim[1]*ff->topGridDim[2];
    // L = 1 + floor(lg(N/M)/3)
    ff->maxLevel = N < 8*M ? 1 : 1 + ilogb(N/M)/3;}
  ff->aCut = (double *)malloc((ff->maxLevel + 1)*sizeof(double));
  ff->aCut[ff->maxLevel] = aL;
  for (int l = ff->maxLevel-1; l >= 0; l--)
    ff->aCut[l] = 0.5*ff->aCut[l+1];
  if (! ff->tolDir){
    if (ff->tolRec) ff->tolDir = ff->tolRec;
    else ff->tolDir = 0.1*pow(0.5*ff->relCutoff, -ff->orderAcc);}
  if (! ff->tolRec)
    ff->tolRec = ff->tolDir;
  
  // build tau(s) = tau_0 + tau_1*s + ... + tau_{ord-1} *s^{ord-1}
  ff->tau = (double *)malloc(ff->orderAcc*sizeof(double));
  ff->tau[0] = 1.;
  for (int i = 1; i < ff->orderAcc; i++)
    ff->tau[i] = - (1. - 0.5/(double)i)*ff->tau[i-1];
  
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

  // build two-scale stencil J[n], n = -nu/2, ..., nu/2
  ff->J = (double *)calloc(nu + 1, sizeof(double));
  // calculate Pascal's triangle
  double *J0 = ff->J + nu/2;
  J0[0] = 1.;
  for (int i = 1; i <= nu/2; i++){
    for (int j = -i; j < i; j++) J0[j] += J0[j+1];
    for (int j = i; j > -i; j--) J0[j] += J0[j-1];}
  for (int i = -nu/2; i <= nu/2; i++) J0[i] *= pow(2., 1 - nu);
  
  // set beta and kmax
  double pi = 4.*atan(1.);
  double h0 = pow(detA/(double)N, 1./3.);
  double const_ = ff->tolDir*aL/h0;
  double beta = 0.; // beta = beta*aL until after iteration
  double res = erfc(beta) - const_;
  double oldres = 2.*res;
  while (fabs(res) < fabs(oldres)){  // Newton-Raphson
    double dres = -2./sqrt(pi)*exp(-beta*beta);
    beta -= res/dres;
    oldres = res;
    res = erfc(beta) - const_;}
  beta = beta/aL;
  ff->beta = beta;
  const_ = sqrt(pi)*ff->tolRec/(2.*beta*h0);
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

  ff->kLim[0] = (int)(kmax/asx);
  ff->kLim[1] = (int)(kmax/asy);
  ff->kLim[2] = (int)(kmax/asz);
  
  // d^i sigma (1-) = d^i gama(1-) + i * d^(i-1) gama(1-)
  // Implementation of A.1 from periodic.pdf (version 20170811)
  double *dsigma = malloc(sizeof(double) * 2*nu);
  double **dgama = (double **)malloc(sizeof(double *)*2*nu);
  for (int i=0 ; i<2*nu; i++)
    dgama[i] = (double *)malloc(sizeof(double)*nu);
  
  for (int k=nu-1;k>=0 ;k--) dgama[0][k] = 1;
  for (int i = 1; i < 2*nu ; i++) {
    dgama[i][nu-1] = 0;
    for (int k = nu - 1 ; k >= 1; k--) {
      dgama[i][k-1] = ((.5 - k)/k) * (2.*i*dgama[i-1][k] +
                                      (i > 1 ? (i*(i-1.0)*dgama[i-2][k]) : 0.0));
    }
  }
  dsigma[0] = dgama[0][0];
  for (int i = 1 ; i <= 2*nu - 1 ; i++)
    dsigma[i] = dgama[i][0] + i * dgama[i-1][0];
  ff->dsigma = dsigma;
  for (int i=0 ; i<2*nu; i++) free(dgama[i]);
  free(dgama);

  if (ff->altSplitting) {
    Matrix A =  *(Matrix *)edges;
    double ax = sqrt(A.xx*A.xx + A.xy*A.xy + A.xz*A.xz);
    double ay = sqrt(A.yx*A.yx + A.yy*A.yy + A.yz*A.yz);
    double az = sqrt(A.zx*A.zx + A.zy*A.zy + A.zz*A.zz);
    double hL = fmax(ax/ff->topGridDim[0],
                fmax(ay/ff->topGridDim[1],az/ff->topGridDim[2]));
    int L = ff->maxLevel ;
    double kmax = 1.0/(hL* pow(nu * pow(2.0,L+nu+1./nu)*ff->tolRec,1./(nu-1.))) ;
    ff->kLim[0] = (int)(kmax/asx);
    ff->kLim[1] = (int)(kmax/asy);
    ff->kLim[2] = (int)(kmax/asz);
  }

  // build anti-blurring operator
  omegap(ff);
  
  // build grid-to-grid stencil
  // assume DFT used
  // allocate memory for khat
  // special treatment for level l = maxLevel
  int l = ff->maxLevel;
  ff->khat = (double **)malloc((l+1)*sizeof(double *));
  int di = ff->topGridDim[0], dj = ff->topGridDim[1], dk = ff->topGridDim[2];
  ff->khat[l] = (double *)calloc(di*dj*dk, sizeof(double));
#ifdef NO_FFT
  ff->fftw_in = (double complex *)malloc(di*dj*dk*sizeof(double complex));
#else
  ff->fftw_in = (fftw_complex *)fftw_malloc(di*dj*dk*sizeof(fftw_complex));
  ff->forward = fftw_plan_dft_3d(di, dj, dk, ff->fftw_in, ff->fftw_in,
                                FFTW_FORWARD, FFTW_ESTIMATE);
  ff->backward = fftw_plan_dft_3d(di, dj, dk, ff->fftw_in, ff->fftw_in,
                                FFTW_BACKWARD, FFTW_ESTIMATE);
#endif

  int dmax = 2*ff->nLim + 1;
  for (l = ff->maxLevel-1; l > 0; l--){
    di *= 2; dj *= 2; dk *= 2;
    int sdi = dmax < di ? dmax : di;
    int sdj = dmax < dj ? dmax : dj;
    int sdk = dmax < dk ? dmax : dk;
    ff->khat[l] = (double *)calloc(sdi*sdj*sdk, sizeof(double));}
  // following is for purpose of debugging
  if (strcmp(o.test, "nobuild") == 0) return;
  // compute stencils ff->khat[l]
  
 
  
  FF_rebuild(ff, edges);
}

// for each computational parameter
// also have a get method FF_set_... that returns the parameter
double FF_get_relCutoff(FF*ff) {
  return ff->relCutoff;}
int FF_get_orderAcc(FF*ff) {
  return ff->orderAcc;}
int FF_get_maxLevel(FF *ff) {
  return ff->maxLevel;}
void FF_get_topGridDim(FF*ff, int topGridDim[3]) {
  *(Triple *)topGridDim = *(Triple *)ff->topGridDim;}
double FF_get_tolDir(FF*ff) {
  return ff->tolDir;}
double FF_get_tolRec(FF*ff) {
  return ff->tolRec;}
bool FF_get_FFT(FF *ff){
        return ff->FFT;}
double FF_get_cutoff(FF *ff) {
  return ff->aCut[0];}
bool FF_get_altSplitting(FF *ff) {
  return ff->altSplitting;}

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
  double a0 = ff->aCut[0];
  return Q2*C*pow(hmin, nu-2)/(pow(detA, 1./3.)*sqrt((double)N)*pow(a0, nu-1));
}


// The rebuild method is called every time the periodic cell changes.
// It initializes edges and calculate the grid2grid stencils.
// helper functions:
static void dALp1(FF *ff, Triple gd, Triple sd, double kh[], double detA);
static void dAL(FF *ff, Triple gd, Triple sd, double kh[], double detA);
static void kaphatA(FF *ff, int l, Triple gd, Triple sd, double kh[],
                    Vector as);
static void DFT(Triple gd, double dL[], double khatL[]);
static void FFT(FF *ff, Triple gd, double dL[], double khatL[]);
void FF_rebuild(FF *ff, double edges[3][3]) {
  *(Matrix *)ff->A = *(Matrix *)edges;
  Matrix A = *(Matrix *)ff->A;
  // calculate coeffs for const part of energy
  double *tau = ff->tau;
  int nu = ff->orderAcc;
  double a_0 = ff->aCut[0];
  int L = ff->maxLevel;
  double aL=ff->aCut[L];
  double beta = ff->beta;
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
        if (x*x + y*y + z*z == 0) continue;
        double rx = A.xx*x + A.xy*y + A.xz*z;
        double ry = A.yx*x + A.yy*y + A.yz*z;
        double rz = A.zx*x + A.zy*y + A.zz*z;
        double r2 = rx*rx + ry*ry + rz*rz;
        if (r2 >= a_0*a_0) continue;
        if (r2 != 0.) ff->coeff1 += 1./sqrt(r2);
        double s = r2/(a_0*a_0) - 1.;
        // gamma(r/a_0) = tau(r2/a_0^2 - 1)
        double gam = tau[nu-1];
        for (int k = nu - 2; k >= 0; k--)
          gam = tau[k] + s*gam;
        ff->coeff1 -= gam/a_0;}
  // compute f->coeff2
  double pi = 4.*atan(1.);
  if (ff->altSplitting)
    ff->coeff2 = pi * aL * aL / ((4 * nu + 2) * detA );
  else
    ff->coeff2 = pi/(beta*beta*detA);  // coeff of 1/2(sum_i q_i)^2
  

  // build grid-to-grid stencil for levels L, L + 1
  // ff->khat[L] = d^{L+1} + DFT of khat^L
  Triple gd = *(Triple *)ff->topGridDim;
  // determine range of kernel evaluations khat^L
  int kdmax = 2*ff->nLim + 1;
  Triple sd = gd;
  sd.x = kdmax < gd.x ? kdmax : gd.x;
  sd.y = kdmax < gd.y ? kdmax : gd.y;
  sd.z = kdmax < gd.z ? kdmax : gd.z;
  // calculate level L kappa hat
  double *kh = ff->khat[L];
  if (ff->altSplitting) {
    dAL(ff, gd, gd, kh, detA);
  } else {
    double *khL = (double *)calloc(sd.x*sd.y*sd.z, sizeof(double));
    kaphatA(ff, L, gd, sd, khL, as); // add in real space contribution
    double *khatL = (double *)calloc(gd.x*gd.y*gd.z, sizeof(double));
    // expand kh into khatL
    for (int mx = - sd.x/2; mx <= (sd.x - 1)/2; mx++)
      for (int my =  - sd.y/2; my <= (sd.y - 1)/2; my++)
        for (int mz =  - sd.z/2; mz <= (sd.z - 1)/2; mz++){
          int nx = (mx + gd.x) % gd.x;
          int ny = (my + gd.y) % gd.y;
          int nz = (mz + gd.z) % gd.z;
          khatL[(nx*gd.y + ny)*gd.z + nz]
            = khL[((mx + sd.x/2)*sd.y + my + sd.y/2)*sd.z + mz + sd.z/2];}
    free(khL);
    if (ff->FFT) FFT(ff, gd, kh, khatL);
    else DFT(gd, kh, khatL);
    free(khatL);
    // add in d^{L+1}(A)
    // d^{L+1}(A)_n = sum_k chi(k) c'(k) exp(2 pi i k . H_L n)
    dALp1(ff, gd, gd, kh, detA);
  }
  // build grid-to-grid stencil for levels L-1, ..., 1
  for (int l = L - 1; l > 0; l--){
    gd.x *= 2; gd.y *= 2; gd.z *= 2;
    sd.x = kdmax < gd.x ? kdmax : gd.x;
    sd.y = kdmax < gd.y ? kdmax : gd.y;
    sd.z = kdmax < gd.z ? kdmax : gd.z;
    double *kh = ff->khat[l];
    for (int i = 0; i < sd.x*sd.y*sd.z; i++) kh[i] = 0.;
    //:::o.time = clock();
    kaphatA(ff, l, gd, sd, kh, as);
    //:::clock_t end = clock();
    //:::printf("l = %d, elapsed time = %f\n",
    //:::      l, (double)(end - o.time)/CLOCKS_PER_SEC);
  }
}

void FF_delete(FF *ff) {
  free(ff->aCut);
  free(ff->tau);
  free(ff->Q);
  free(ff->J);
  for (int l = 1; l <= ff->maxLevel; l++){
    free(ff->khat[l]);
    for (int alpha = 0; alpha < 3; alpha++)
      free(ff->omegap[l][alpha]);}
  free(ff->khat);
  free(ff->omegap);
  free(ff->dsigma);
  for (int alpha = 0; alpha < 3; alpha++)
    free(ff->cL[alpha]);
#ifdef NO_FFT
  free(ff->fftw_in);
#else
  fftw_destroy_plan(ff->forward);
  fftw_destroy_plan(ff->backward);
  fftw_free(ff->fftw_in);
#endif
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
  // construct ff->omegap and ff->cL
  int nu = ff->orderAcc;
  // nu positive even integer
  int nuby2 = nu/2;

  // Compute Phi(n), n = 0, 1,..., nuby2-1
  double *Phi = (double *)calloc(nuby2, sizeof(double));
  for(int i = 1; i < nuby2; i++)  Phi[nuby2-i] = ff->Q[i*nu];
  for (int k = 0; k < nu; k++) Phi[0] += ff->Q[(nuby2-1)*nu + k];

  int M[3] = {ff->topGridDim[0], ff->topGridDim[1], ff->topGridDim[2]};
  double pi = 4.*atan(1.);
  ff->omegap = (double *(*)[3])malloc((ff->maxLevel+1)*sizeof(double *[3]));
  int opLim = 2*ff->nLim;
  for (int l = ff->maxLevel; l >= 1; l--){
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
      int opD = Ma/2 + 1;
      if (Ma/2 > opLim) opD = opLim + 1;
      ff->omegap[l][alpha] = (double *)calloc(opD, sizeof(double));
      for (int m = 0; m < opD; m++){
        double om_m = c[0];
        for (int k = 1; k <= (Ma - 1)/2; k++)
          om_m += 2.*c[k]*cos(2.*pi*(double)(k*m)/(double)Ma);
        if (Ma%2 == 0) om_m += c[Ma/2]*cos(pi*(double)(m));
        ff->omegap[l][alpha][m] = om_m/(double)Ma;}
      if (l == ff->maxLevel) ff->cL[alpha] = c;
      else free(c);
    }
    M[0] *= 2, M[1] *= 2, M[2] *= 2;}
  free(Phi);
}

static void dAL(FF *ff, Triple gd, Triple sd, double kh[], double detA){
  // add d^{L+1}(A) to kh
  Matrix Ai = *(Matrix *)ff-> Ai;
  int nu = ff->orderAcc;
  double *sigmad = ff->dsigma;
  int L = ff->maxLevel;
  double aL = ff->aCut[L];
  double pi = 4.*atan(1.);
  // loop on vec k
  // for kx = 0, 1, -1, ..., kLim.x, -kLim.x
  int klimx, klimy, klimz;
  if (ff->kLimUserSpecified >= 0)
    klimx = klimy = klimz = ff->kLimUserSpecified ;
  else
    klimx = ff->kLim[0], klimy = ff->kLim[1], klimz = ff->kLim[2];
  for (int kx = - klimx ; kx <= klimx; kx++){
    int kx1 = (kx % gd.x + gd.x) % gd.x;
    int kx0 = (kx1 + gd.x/2) % gd.x - gd.x/2;
    double cLx = ff->cL[0][abs(kx0)];
    for (int ky = - klimy; ky <= klimy; ky++){
      int ky1 = (ky % gd.y + gd.y) % gd.y;
      int ky0 = (ky1 + gd.y/2) % gd.y - gd.y/2;
      double cLxy = cLx*ff->cL[1][abs(ky0)];
      for (int kz = - klimz; kz <= klimz; kz++){
        int kz1 = (kz % gd.z + gd.z) % gd.z;
        int kz0 = (kz1 + gd.z/2) % gd.z - gd.z/2;
        double cLxyz = cLxy*ff->cL[2][abs(kz0)];
        double fkx = (double)kx, fky = (double)ky, fkz = (double)kz;
        double Aikx = Ai.xx*fkx + Ai.yx*fky + Ai.zx*fkz,
        Aiky = Ai.xy*fkx + Ai.yy*fky + Ai.zy*fkz,
        Aikz = Ai.xz*fkx + Ai.yz*fky + Ai.zz*fkz;
        double k2 = Aikx*Aikx + Aiky*Aiky + Aikz*Aikz;
        double chi_cL = 0.0;
        if (k2 != 0) {
          double k = sqrt(k2);
          double chisum = 0;
          for (int j = nu/2 ; j<=nu -1 ;j++)
            chisum += pow(-1,j) *
              (-cos(pi*k*aL) * sigmad[2*j]   / pow(pi*k*aL,2*j+1)
               +sin(pi*k*aL) * sigmad[2*j+1] / pow(pi*k*aL,2*j+2));
          chisum *= aL/(k*detA) ;
          chi_cL = chisum * cLxyz;
        }
        kh[(kx1*gd.y + ky1)*gd.z + kz1]
        += chi_cL*gd.x*gd.y*gd.z;
      }
    }
  }
}

static void dALp1(FF *ff, Triple gd, Triple sd, double kh[], double detA){
  // add d^{L+1}(A) to kh
  Matrix Ai = *(Matrix *)ff-> Ai;
  double pi = 4.*atan(1.);
  double pidetA = pi*fabs(detA);
  double pi2beta2 = pow(pi/ff->beta, 2);
  // loop on vec k
  // for kx = 0, 1, -1, ..., kLim.x, -kLim.x
  int klimx, klimy, klimz;
  if (ff->kLimUserSpecified >= 0)
    klimx = klimy = klimz = ff->kLimUserSpecified ;
  else
    klimx = ff->kLim[0], klimy = ff->kLim[1], klimz = ff->kLim[2];
  for (int kx = - klimx ; kx <= klimx; kx++){
    int kx1 = (kx % gd.x + gd.x) % gd.x;
    int kx0 = (kx1 + gd.x/2) % gd.x - gd.x/2;
    double cLx = ff->cL[0][abs(kx0)];
    for (int ky = - klimy; ky <= klimy; ky++){
      int ky1 = (ky % gd.y + gd.y) % gd.y;
      int ky0 = (ky1 + gd.y/2) % gd.y - gd.y/2;
      double cLxy = cLx*ff->cL[1][abs(ky0)];
      for (int kz = - klimz; kz <= klimz; kz++){
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

static void DFT(Triple gd, double dL[], double khatL[]){
  double twopi = 8.*atan(1.);
  for (int kx = 0; kx < gd.x; kx++)
    for (int ky = 0; ky < gd.y; ky++)
      for (int kz = 0; kz < gd.z; kz++){
        double dLk = 0.;
        double cx = twopi*(double)kx/gd.x;
        double cy = twopi*(double)ky/gd.y;
        double cz = twopi*(double)kz/gd.z;
        for (int nx = 0; nx < gd.x; nx++)
          for (int ny = 0.; ny < gd.y; ny++)
            for (int nz = 0.; nz < gd.z; nz++)
              dLk += cos(cx*(double)nx + cy*(double)ny + cz*(double)nz)
                *khatL[(nx*gd.y + ny)*gd.z + nz];
        dL[(kx*gd.y + ky)*gd.z + kz] = dLk;
      }
}

static void FFT(FF *ff, Triple gd, double dL[], double khatL[]){
#ifdef NO_FFT
  ;
#else
  // dL = DFT of khatL
  for (int i = 0; i < gd.x*gd.y*gd.z; i++)
    ff->fftw_in[i] = (fftw_complex)khatL[i];
  fftw_execute(ff->forward);
  for (int i = 0; i < gd.x*gd.y*gd.z; i++)
    dL[i] = creal(ff->fftw_in[i]);
#endif
}

static double kappaA(FF *ff, int l, Vector s, Vector as);  // kappa_l(A s; A)
static void kaphatA(FF *ff, int l, Triple gd, Triple sd, double kh[],
                    Vector as){
  // add kappaHat_l to kh
  //clock_t begin = clock();
  int kdmax = 2*ff->nLim + 1;
  int kdx = kdmax < gd.x ? kdmax : gd.x;
  int kdy = kdmax < gd.y ? kdmax : gd.y;
  int kdz = kdmax < gd.z ? kdmax : gd.z;
  // construct array of kappa values
  double *kap = (double *)malloc(kdx*kdy*kdz*sizeof(double));
  Vector s;
  for (int i = - kdx/2; i <= (kdx - 1)/2; i++){
    double *kapi = kap + (i + kdx/2)*kdy*kdz;
    s.x = (double)i/(double)gd.x;
    for (int j = - kdy/2; j <= (kdy - 1)/2; j++){
      double *kapij = kapi + (j + kdy/2)*kdz;
      s.y = (double)j/(double)gd.y;
      for (int k = - kdz/2; k <= (kdz - 1)/2; k++){
        s.z = (double)k/(double)gd.z;
        kapij[k + kdz/2] = kappaA(ff, l, s, as);}}}
  if (strcmp(o.test, "test_kappaA") == 0)
    for (int i = 0; i < kdx*kdy*kdz; i++) o.kappa[l][i] = kap[i];
  double **op = ff->omegap[l];
  int opLim = 2*ff->nLim;
  //clock_t end = clock();
  //-printf("elapsed time = %f, iterations = %d\n",
  //-(double)(end - o.time)/CLOCKS_PER_SEC, kdx*kdy*kdz);
  //begin = clock();
  //construct kappa hat element by element
  if (l == ff->maxLevel)
    for (int i1 = - sd.x/2; i1 <= (sd.x - 1)/2; i1++){
      double *khi = kh + (i1 + sd.x/2)*sd.y*sd.z;
      for (int j1 = - sd.y/2; j1 <= (sd.y - 1)/2; j1++){
        double *khij = khi + (j1 + sd.y/2)*sd.z;
        for (int k1 = - sd.z/2; k1 <= (sd.z - 1)/2; k1++){
          double khijk = 0.;
          for (int i0 = - kdx/2; i0 <= (kdx - 1)/2; i0++){
            int i = (i0 - i1 + (3*gd.x)/2)%gd.x - gd.x/2;
            double opi = op[0][abs(i)];
            for (int j0 = - kdy/2; j0 <= (kdy - 1)/2; j0++){
              int j = (j0 - j1 + (3*gd.y)/2)%gd.y - gd.y/2;
              double opij = opi*op[1][abs(j)];
              for (int k0 = - kdz/2; k0 <= (kdz - 1)/2; k0++){
                int k = (k0 - k1 + (3*gd.z)/2)%gd.z - gd.z/2;
                double opijk = opij*op[2][abs(k)];
                double kapijk
                  = kap[((i0 + kdx/2)*kdy + j0 + kdy/2)*kdz + k0 + kdz/2];
                khijk += opijk*kapijk;
              }}}
          khij[k1 + sd.z/2] = khijk;}}}
  else
    for (int i1 = - sd.x/2; i1 <= (sd.x - 1)/2; i1++){
      double *khi = kh + (i1 + sd.x/2)*sd.y*sd.z;
      for (int j1 = - sd.y/2; j1 <= (sd.y - 1)/2; j1++){
        double *khij = khi + (j1 + sd.y/2)*sd.z;
        for (int k1 = - sd.z/2; k1 <= (sd.z - 1)/2; k1++){
          double khijk = 0.;
          for (int i0 = - kdx/2; i0 <= (kdx - 1)/2; i0++){
            int i = (i0 - i1 + (3*gd.x)/2)%gd.x - gd.x/2;
            double opi = op[0][abs(i)];
            for (int j0 = - kdy/2; j0 <= (kdy - 1)/2; j0++){
              int j = (j0 - j1 + (3*gd.y)/2)%gd.y - gd.y/2;
              double opij = opi*op[1][abs(j)];
              for (int k0 = - kdz/2; k0 <= (kdz - 1)/2; k0++){
                int k = (k0 - k1 + (3*gd.z)/2)%gd.z - gd.z/2;
                double opijk = opij*op[2][abs(k)];
                double kapijk
                  = kap[((i0 + kdx/2)*kdy + j0 + kdy/2)*kdz + k0 + kdz/2];
                khijk += opijk*kapijk;
              }}}
          khij[k1 + sd.z/2] = khijk;}}}
  free(kap);
  //end = clock();
  //-printf("elapsed time = %f, iterations = %d\n",
  //-(double)(end - o.time)/CLOCKS_PER_SEC, sd.x*sd.y*sd.z);
}

static double kappaA(FF *ff, int l, Vector s, Vector as){
  // kappa_l(A s; A)  ::: might separate l < L and l = L
  Matrix A = *(Matrix *)ff->A;
  double *tau = ff->tau;
  int nu = ff->orderAcc;
  double a_l = ff->aCut[l];
  double beta = ff->beta;
  double rootPi = sqrt(4.*atan(1.));

  double pxmin = floor(s.x - a_l*as.x + 1.), pxmax = ceil(s.x + a_l*as.x - 1.);
  double pymin = floor(s.y - a_l*as.y + 1.), pymax = ceil(s.y + a_l*as.y - 1.);
  double pzmin = floor(s.z - a_l*as.z + 1.), pzmax = ceil(s.z + a_l*as.z - 1.);
  double kap_s = 0.;
  for (double px = pxmin; px <= pxmax; px++)
    for (double py = pymin; py <= pymax; py++)
      for (double pz = pzmin; pz <= pzmax; pz++){ // could do just pz >= 0
        double x = s.x - px, y = s.y - py, z = s.z - pz;
        double rx = A.xx*x + A.xy*y + A.xz*z;
        double ry = A.yx*x + A.yy*y + A.yz*z;
        double rz = A.zx*x + A.zy*y + A.zz*z;
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        double rho = r/a_l;
        if (rho >= 1.) continue;
        // one could instead precompute a list of Vectors (p_i, p_j, p_k)
        double gam, gam2;
        if (l < ff->maxLevel){
          if (rho >= 1.) gam = 1./rho;
          else{ // gam = tau(rho^2 - 1)
            double s = rho*rho - 1.;
            gam = tau[nu-1];
            for (int k = nu-2; k >= 0; k--)
              gam = tau[k] + s*gam;}
          gam /= a_l;}
        else
          gam = r != 0. ? erf(beta*r)/r : 2.*beta/rootPi;
        rho *= 2.;
        if (rho >= 1.) gam2 = 1./rho;
        else{
          double s = rho*rho - 1.;
          gam2 = tau[nu-1];
          for (int k = nu-2; k >= 0; k--)
            gam2 = tau[k] + s*gam2;}
        gam2 *= 2./a_l;
        kap_s += gam2 - gam;}
  return kap_s;}

//::: not used
double kappa(FF *ff, int l, double s[3], double edges[3][3]){
  Matrix Ai = *(Matrix *)edges;
  double detA = invert(&Ai);
  *(Matrix *)ff->Ai = Ai;
  Vector as = {sqrt(Ai.xx*Ai.xx + Ai.yx*Ai.yx + Ai.zx*Ai.zx),
               sqrt(Ai.xy*Ai.xy + Ai.yy*Ai.yy + Ai.zy*Ai.zy),
               sqrt(Ai.xz*Ai.xz + Ai.yz*Ai.yz + Ai.zz*Ai.zz)};
  return kappaA(ff, l, *(Vector *)s, as);}
