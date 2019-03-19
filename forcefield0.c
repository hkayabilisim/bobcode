// file forcefield.c
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
#include "forcefield.h"

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
    double asprod = asx * asy * asz;
    double Ms;
    if (ff->maxLevel)
      Ms = pow(2.,(double)ff->maxLevel)*(double)N;
    else {
      Ms = sqrt((double)N);
      double asmin = fmin(asx, fmin(asy, asz));
      double Msmin = pow(0.5*ff->relCutoff*asmin, 3)/asprod;
      Ms = fmax(Ms, Msmin);}
    double numerator = pow(Ms*asprod, 1./3.);
    ff->topGridDim[0] = (int)ceil(numerator/asx);
    ff->topGridDim[1] = (int)ceil(numerator/asy);
    ff->topGridDim[2] = (int)ceil(numerator/asz);}
  double aL = ff->relCutoff/fmax(ff->topGridDim[0]*asx,
              fmax(ff->topGridDim[1]*asy,ff->topGridDim[2]*asz));
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
  int dmax = 2*ff->nLim + 1;
  for (l = ff->maxLevel-1; l > 0; l--){
    di *= 2; dj *= 2; dk *= 2;
    int sdi = dmax < di ? dmax : di;
    int sdj = dmax < dj ? dmax : dj;
    int sdk = dmax < dk ? dmax : dk;
    sdk = sdk/2 + 1; // use the fact that khat_{i,j,k} = khat_{i,j,-k}
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
static void kaphatALp1(FF *ff, Triple gd, Triple sd, double kh[], double detA);
static void dALp1(FF *ff, Triple gd, Triple sd, double kh[], double detA);
static void kaphatA(FF *ff, int l, Triple gd, Triple sd, double kh[],
                    Vector as);
static void	DFT(Triple gd, double dL[], double khatL[]);
void FF_rebuild(FF *ff, double edges[3][3]) {
  *(Matrix *)ff->A = *(Matrix *)edges;
  Matrix A = *(Matrix *)ff->A;
  // calculate coeffs for const part of energy
  double *tau = ff->tau;
  int nu = ff->orderAcc;
  double a_0 = ff->aCut[0];
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
  ff->coeff2 = pi/(beta*beta*detA);  // coeff of 1/2(sum_i q_i)^2

  // build grid-to-grid stencil for level L
  Triple gd = *(Triple *)ff->topGridDim;
  // determine range of kernel evaluations
  Triple sd;
  sd.x = gd.x, sd.y = gd.y, sd.z = gd.z/2 + 1;  // stencil dimensions
   // using the fact that khat_{i,j,k} = khat_{i,j,-k}
  // calculate level L+1 kappa hat
  int L = ff->maxLevel;
  double *kh = ff->khat[L];
  // kappa_hat_n = sum_k chi(k) c'(k) exp(2 pi i k . H_L n)
  // kh = d^{L+1} + DFT of khat^L
	for (int i = 0; i < sd.x*sd.y*sd.z; i++) kh[i] = 0.;
	kaphatA(ff, L, gd, sd, kh, as); // add in real space contribution
	double *khatL = (double *)malloc(gd.x*gd.y*gd.z*sizeof(double));
	// expand kh into khatL
	for (int i = 0; i < gd.x; i++)
		for (int j = 0; j < gd.y; j++)
			for (int k = 0; k < gd.z; k++){
				int i_, j_, k_;
				if (k < (gd.z + 1)/2)
					i_ = gd.x - i, j_ = gd.y - j, k_ = - k;
				else
					i_ = i, j_ = j, k_ = k - gd.z;
				i_ = (i_ + gd.x/2)%gd.x - gd.x/2;
				j_ = (j_ + gd.y/2)%gd.y - gd.y/2;
				khatL[(i*gd.y + j)*gd.z + k]
					= kh[((i_ + sd.x/2)*sd.y + j_ + sd.y/2)*sd.z + k_ + sd.z - 1];
			}
	double *dL = (double *)malloc(gd.x*gd.y*gd.z*sizeof(double));
	DFT(gd, dL, khatL);
	free(khatL);
	// compress dL into kh
	for (int i_ = - gd.x/2; i_ <= (gd.x - 1)/2; i_++)
		for (int j_ = - gd.y/2; j_ <= (gd.y - 1)/2; j_++)
			for (int k_ = - gd.z/2; k_ <= 0; k_++){
				int i = (i_ + gd.x)%gd.x, j	= (j_ + gd.y)%gd.y, k	= (k_ + gd.z)%gd.z;
				kh[((i_ + sd.x/2)*sd.y + j_ + sd.y/2)*sd.z + k_ + sd.z - 1]
					= dL[(i*gd.y + j)*gd.z + k];
			}
	free(dL);
	dALp1(ff, gd, sd, kh, detA);  // add in d^{L+1}(A)
  
  // build grid-to-grid stencil for levels L-1, ..., 1
  int kdmax = 2*ff->nLim + 1;
  for (int l = L - 1; l > 0; l--){
    double *kh = ff->khat[l];
    for (int i = 0; i < sd.x*sd.y*sd.z; i++) kh[i] = 0.;
    gd.x *= 2; gd.y *= 2; gd.z *= 2;
    sd.x = kdmax < gd.x ? kdmax : gd.x;
    sd.y = kdmax < gd.y ? kdmax : gd.y;
    sd.z = kdmax < gd.z ? kdmax : gd.z;
    sd.z = sd.z/2 + 1; // using the fact that khat_{i,j,k} = khat_{i,j,-k}
    kaphatA(ff, l, gd, sd, kh, as);
  }
}

void FF_delete(FF *ff) {
  ;
  //... free memory;
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

static void kaphatALp1(FF *ff, Triple gd, Triple sd, double kh[], double detA){
	// set kh to kappaHat^{L+1}(A)
  // determine range of kernel evaluations
  Matrix Ai = *(Matrix *)ff-> Ai;
  Triple kLim = *(Triple *)ff->kLim;
  // loop on vec k
  double pi = 4.*atan(1.);
  double pidetA = pi*fabs(detA);
  double pi2beta2 = pow(pi/ff->beta, 2);
  // for kx = 0, 1, -1, ..., kLim.x, -kLim.x
  int dx = 2*kLim.x + 1, dy = 2*kLim.y + 1, dz = 2*kLim.z + 1;
  double *d = (double *)calloc(dx*dy*dz, sizeof(double));
  for (int kx = - kLim.x ; kx <= kLim.x; kx++){
    int kx0 = ((kx + gd.x/2) % gd.x + gd.x) % gd.x - gd.x/2;
    double cLx = ff->cL[0][abs(kx0)];
    for (int ky = - kLim.y; ky <= kLim.y; ky++){
      int ky0 = ((ky + gd.y/2) % gd.y + gd.y) % gd.y - gd.y/2;
      double cLxy = cLx*ff->cL[1][abs(ky0)];
      for (int kz = - kLim.z; kz <= kLim.z; kz++){
        int kz0 = ((kz + gd.z/2) % gd.z + gd.z) % gd.z - gd.z/2;
        double cLxyz = cLxy*ff->cL[2][abs(kz0)];
        double fkx = (double)kx, fky = (double)ky, fkz = (double)kz;
        double Aikx = Ai.xx*fkx + Ai.yx*fky + Ai.zx*fkz,
          Aiky = Ai.xy*fkx + Ai.yy*fky + Ai.zy*fkz,
          Aikz = Ai.xz*fkx + Ai.yz*fky + Ai.zz*fkz;
        double k2 = Aikx*Aikx + Aiky*Aiky + Aikz*Aikz;
        double chi_cL = k2 == 0 ? 0 : exp(- pi2beta2*k2)/(pidetA*k2)*cLxyz;
        d[((kx0 + dx/2)*dy + ky0 + dy/2)*dz + kz0 + dz/2] += chi_cL;
      }
    }
  }
  
  for (int i = 0; i < sd.x*sd.y*sd.z; i++) kh[i] = 0.;
  double complex zetax = cexp(2.*pi*I/(double)gd.x);
  // for kx = 0, 1, -1, ..., kLim.x, -kLim.x
  double complex zetax_k = 1.;
  for (int kx = 0 ; kx != kLim.x + 1
         ; zetax_k =  kx > 0 ? conj(zetax_k) : zetax*conj(zetax_k),
         kx = kx > 0 ? - kx : 1 - kx){
    int kx0 = abs(((kx + gd.x/2) % gd.x + gd.x) % gd.x - gd.x/2);
    double cLx = ff->cL[0][kx0];
    double complex zetay = cexp(2.*pi*I/(double)gd.y);
    double complex zetay_k = 1.;
    for (int ky = 0; ky != kLim.y + 1
           ; zetay_k =  ky > 0 ? conj(zetay_k) : zetay*conj(zetay_k),
           ky = ky > 0 ? - ky : 1 - ky){
      int ky0 = abs(((ky + gd.y/2) % gd.y + gd.y) % gd.y - gd.y/2);
      double cLxy = cLx*ff->cL[1][ky0];
      double complex zetaz = cexp(2.*pi*I/(double)gd.z);
      double complex zetaz_k = 1.;
      for (int kz = 0; kz != kLim.z + 1
             ; zetaz_k =  kz > 0 ? conj(zetaz_k) : zetaz*conj(zetaz_k),
             kz = kz > 0 ? - kz : 1 - kz){
        int kz0 = abs(((kz + gd.z/2) % gd.z + gd.z) % gd.z - gd.z/2);
        double cLxyz = cLxy*ff->cL[2][kz0];
        double chi_cL = d[((kx + dx/2)*dy + ky + dy/2)*dz + kz + dy/2];
        // distribute chi(k; A)c(k)^2
        // for - sd.x/2 <= nx <= (sd.x - 1)/2
        double complex zetax_kn = 1.;
        for (int nx = 0; nx != sd.x/2 + 1
               ; zetax_kn = nx > 0 ? conj(zetax_kn) : zetax_k*conj(zetax_kn),
               nx = nx > 0 ? - nx : 1 - nx){
          if (nx == (sd.x - 1)/2 + 1) continue;
          double *khi = kh + (nx+sd.x/2)*sd.y*sd.z;
          // for - sd.y/2 <= ny <= (sd.y - 1)/2
          double complex zetay_kn = 1.;
          for (int ny = 0; ny != sd.y/2 + 1
                 ; zetay_kn = ny > 0 ? conj(zetay_kn) : zetay_k*conj(zetay_kn),
                 ny = ny > 0 ? - ny : 1 - ny){
            if (ny == (sd.y - 1)/2 + 1) continue;
            double *khij = khi + (ny+sd.y/2)*sd.z;
            // for 1 - sd.z <= nz <= 0
            double complex zetaz_kn = 1.;
            for (int nz = 0; nz != - sd.z
                   ; zetaz_kn *= conj(zetaz_k), nz--){
              // add chi_cL*zetax_kn*zetay_kn*zetaz_kn to kh_ijk
              khij[nz + sd.z - 1] += chi_cL*creal(zetax_kn*zetay_kn*zetaz_kn);
              if (strcmp(o.test, "test_kaphatA") == 0){
                int offset = (int)(khij + nz + sd.z - 1 - kh);
                o.khatLp1[offset] += chi_cL*creal(zetax_kn*zetay_kn*zetaz_kn);}
              if (strcmp(o.test, "test_kappaA") == 0){
                int offset = (int)(khij + nz + sd.z - 1 - kh);
                int L = ff->maxLevel;
                o.kappa[L+1][offset]
                  += chi_cL*creal(zetax_kn*zetay_kn*zetaz_kn)/cLxyz;}
            }
          }
        }
      }
    }
  }
  free(d);
}

static void dALp1(FF *ff, Triple gd, Triple sd, double kh[], double detA){
  // add d^{L+1}(A) to kh
  Matrix Ai = *(Matrix *)ff-> Ai;
  Triple kLim = *(Triple *)ff->kLim;
  double pi = 4.*atan(1.);
  double pidetA = pi*fabs(detA);
  double pi2beta2 = pow(pi/ff->beta, 2);
  // loop on vec k
  // for kx = 0, 1, -1, ..., kLim.x, -kLim.x
  for (int kx = - kLim.x ; kx <= kLim.x; kx++){
    int kx0 = ((kx + gd.x/2) % gd.x + gd.x) % gd.x - gd.x/2;
    double cLx = ff->cL[0][abs(kx0)];
    for (int ky = - kLim.y; ky <= kLim.y; ky++){
      int ky0 = ((ky + gd.y/2) % gd.y + gd.y) % gd.y - gd.y/2;
      double cLxy = cLx*ff->cL[1][abs(ky0)];
      for (int kz = - kLim.z; kz <= kLim.z; kz++){
        int kz0 = ((kz + gd.z/2) % gd.z + gd.z) % gd.z - gd.z/2;
        double cLxyz = cLxy*ff->cL[2][abs(kz0)];
        double fkx = (double)kx, fky = (double)ky, fkz = (double)kz;
        double Aikx = Ai.xx*fkx + Ai.yx*fky + Ai.zx*fkz,
          Aiky = Ai.xy*fkx + Ai.yy*fky + Ai.zy*fkz,
          Aikz = Ai.xz*fkx + Ai.yz*fky + Ai.zz*fkz;
        double k2 = Aikx*Aikx + Aiky*Aiky + Aikz*Aikz;
        double chi_cL = k2 == 0 ? 0 : exp(- pi2beta2*k2)/(pidetA*k2)*cLxyz;
				if (kz0 <= 0)
					kh[((kx0 + sd.x/2)*sd.y + ky0 + sd.y/2)*sd.z + kz0 + sd.z - 1]
						+= chi_cL*gd.x*gd.y*gd.z;
      }
    }
  }
}

static void	DFT(Triple gd, double dL[], double khatL[]){
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

static double kappaA(FF *ff, int l, Vector s, Vector as);  // kappa_l(A s; A)
static void kaphatA(FF *ff, int l, Triple gd, Triple sd, double kh[],
                    Vector as){
  // add kappaHat_l to kh
  int kdmax = 2*ff->nLim + 1;
  int kdi = kdmax < gd.x ? kdmax : gd.x;
  int kdj = kdmax < gd.y ? kdmax : gd.y;
  int kdk = kdmax < gd.z ? kdmax : gd.z;
  // construct array of kappa values
  double *kap = (double *)malloc(kdi*kdj*kdk*sizeof(double));
  Vector s;
  for (int i = - kdi/2; i <= (kdi - 1)/2; i++){
    double *kapi = kap + (i + kdi/2)*kdj*kdk;
    s.x = (double)i/(double)gd.x;
    for (int j = - kdj/2; j <= (kdj - 1)/2; j++){
      double *kapij = kapi + (j + kdj/2)*kdk;
      s.y = (double)j/(double)gd.y;
      for (int k = - kdk/2; k <= (kdk - 1)/2; k++){
        s.z = (double)k/(double)gd.z;
        kapij[k + kdk/2] = kappaA(ff, l, s, as);}}}
  if (strcmp(o.test, "test_kappaA") == 0)
    for (int i = 0; i < kdi*kdj*kdk; i++) o.kappa[l][i] = kap[i];
  double **op = ff->omegap[l];
  int opLim = 2*ff->nLim;
  //construct kappa hat element by element
  for (int i = - sd.x/2; i <= (sd.x - 1)/2; i++){
    double *khi = kh + (i + sd.x/2)*sd.y*sd.z;
    for (int j = - sd.y/2; j <= (sd.y - 1)/2; j++){
      double *khij = khi + (j + sd.y/2)*sd.z;
      for (int k = 1 - sd.z; k <= 0; k++){
        double khijk = 0.;
        int ioMin = - gd.x/2;
        int ioMax = (gd.x - 1)/2;
        if (gd.x/2 > opLim) ioMin = - opLim, ioMax = opLim;
        for (int io = ioMin; io <= ioMax; io++){
          // calculate i + io modulo gd.x
          int ik = (i + io + gd.x/2 + gd.x)%gd.x - gd.x/2;
          assert( -gd.x/2<= ik && ik<=(gd.x-1)/2 );
          if (ik < - kdi/2 || ik > (kdi-1)/2) continue;
          double *kapi = kap + (ik + kdi/2)*kdj*kdk;
          double opi = op[0][abs(io)];
          int joMin = - gd.y/2;
          int joMax = (gd.y - 1)/2;
          if (gd.y/2 > opLim) joMin = - opLim, joMax = opLim;
          for (int jo = joMin; jo <= joMax; jo++){
            // calculate j + jo modulo gd.y
            int jk = (j + jo + gd.y/2 + gd.y)%gd.y - gd.y/2;
            assert( -gd.y/2<= jk && jk<=(gd.y-1)/2 );
            if (jk < - kdj/2 || jk > (kdj-1)/2) continue;
            double *kapij = kapi + (jk + kdj/2)*kdk;
            double opij = opi*op[1][abs(jo)];
            int koMin = - gd.z/2;
            int koMax = (gd.z - 1)/2;
            if (gd.z/2 > opLim) koMin = - opLim, koMax = opLim;
            for (int ko = koMin; ko <= koMax; ko++){
              // calculate k + ko modulo gd.z
              int kk = (k + ko + gd.z/2 + gd.z)%gd.z - gd.z/2;
              assert( -gd.z/2<= kk && kk<=(gd.z-1)/2 );
              if (kk < - kdk/2 || kk > (kdk-1)/2) continue;
              khijk += opij*op[2][abs(ko)]*kapij[kk + kdk/2];}}}
        khij[k + sd.z - 1] += khijk;}}}
  free(kap);
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
