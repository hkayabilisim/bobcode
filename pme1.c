// file forcefield1.c
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
#include "pme.h"

typedef struct Vector {double x, y, z;} Vector;
typedef struct Matrix {double xx, xy, xz, yx, yy, yz, zx, zy, zz;} Matrix;
typedef struct Triple {int x, y, z;} Triple;

// particles are not assumed to be in periodic cell.
// helper functions
static double partcl2partcl(FF *ff, int N, Vector *force, Vector *position,
                            double *charge);
static void anterpolate(FF *ff, Triple gd, double *q, int N, double *charge,
                        Vector *r);
static void FFTg2g(FF *ff, Triple gd, double *el, double *ql, double *kh);
static void interpolate(FF *ff, int N, Vector *E, Vector *r, Triple gd,
                          double *el);
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

double FF_energy(FF *ff, int N, double (*force)[3], double (*position)[3],
                 double *charge, double *weight) {
  // N may change in case of grand canonical simulations
  Vector *r = (Vector *)position;
  Vector *F = (Vector *)force;
  double *wt = (double *)malloc((1 + 1)*sizeof(double));
  for (int l = 0; l <= 1; l++)
    wt[l] = weight == NULL ? 1. : weight[l];
  // compute forces and energy
  // calculate the constant part of energy
  double qsum2 = 0, q2sum = 0.;
  for (int i = 0; i < N; i++)
    qsum2 += charge[i], q2sum += charge[i]*charge[i];
  qsum2 *= qsum2;
  double energy = 0.5*q2sum*ff->coeff1 - 0.5*qsum2*ff->coeff2;
  // particle-to-particle
  energy *= wt[0];
  msm4g_tic();
  energy += wt[0]*partcl2partcl(ff, N, F, r, charge);
  ff->time_partcl2partcl = msm4g_toc();
  for (int i = 0; i < N; i++)
    F[i].x *= wt[0], F[i].y *= wt[0], F[i].z *= wt[0];
  double **q = (double **)malloc((1 + 1)*sizeof(double *));
  Triple gd = *(Triple *)ff->topGridDim;
  q[1] = (double *)calloc(gd.x*gd.y*gd.z, sizeof(double));
  // calculate q[1]
  anterpolate(ff, gd, q[1], N, charge, r);
  double *el = (double *)calloc(gd.x*gd.y*gd.z, sizeof(double));
  // set e^L = calK^L q^L
  Triple sd = gd;
  for (int m = 0; m < gd.x*gd.y*gd.z; m ++) q[1][m] *= wt[1];
  FFTg2g(ff, gd, el, q[1], ff->khat);
  free(wt);
  // add in grid-level energy
  double grid_en = 0.;
  for (int m = 0; m < gd.x*gd.y*gd.z; m++)
    grid_en += q[1][m]*el[m];
  energy += 0.5*grid_en;
  free(q[1]);
  free(q);
  // initialize grid-level electric field
  Vector *E = (Vector *)malloc(N*sizeof(Vector));
  interpolate(ff, N, E, r, gd, el);
  free(el);
  for (int i = 0; i < N; i++){
    F[i].x += charge[i]*E[i].x;
    F[i].y += charge[i]*E[i].y;
    F[i].z += charge[i]*E[i].z;}
  free(E);
  return energy;
}

static Vector prod(Matrix m, Vector v);
static Vector prodT(Matrix m, Vector v);
static double partcl2partcl(FF *ff, int N, Vector *force, Vector *position,
                            double *charge){
  // return particle-level energy, calculate forces
  Matrix A = *(Matrix *)ff->A;
  Matrix Ai = *(Matrix *)ff->Ai;
  int nu = ff->orderAcc;
  double a_0 = ff->cutoff;
  double energy = 0.;
  for (int i = 0; i < N; i++) force[i].z = force[i].y = force[i].x = 0.;
  Vector as = {sqrt(Ai.xx*Ai.xx + Ai.yx*Ai.yx + Ai.zx*Ai.zx),
               sqrt(Ai.xy*Ai.xy + Ai.yy*Ai.yy + Ai.zy*Ai.zy),
               sqrt(Ai.xz*Ai.xz + Ai.yz*Ai.yz + Ai.zz*Ai.zz)};
  double pxlim = ceil(a_0*as.x - 0.5);
  double pylim = ceil(a_0*as.y - 0.5);
  double pzlim = ceil(a_0*as.z - 0.5);
  // create hash table
  Triple gd = *(Triple *)ff->topGridDim;
  int gd_prod = gd.x*gd.y*gd.z;
  int *first = (int *)malloc(gd_prod*sizeof(int));
  for (int m = 0; m < gd_prod; m++) first[m] = -1;
  int *next = (int *)malloc(N*sizeof(int));
  for (int i = 0; i < N; i++){
    Vector ri = position[i];
    Vector s = prod(Ai, ri);
    s.x = s.x - floor(s.x);
    s.y = s.y - floor(s.y);
    s.z = s.z - floor(s.z);
    Vector t = {(double)gd.x*s.x, (double)gd.y*s.y, (double)gd.z*s.z};
    // Triple m = {(int)floor(t.x), (int)floor(t.y), (int)floor(t.z)};
    Triple m = {(int)t.x, (int)t.y, (int)t.z};
    int m_ = (m.x*gd.y + m.y)*gd.z + m.z;
    next[i] = first[m_];
    first[m_] = i;
  }
  // compute ranges - can be part of preprocessing
  Vector a = {sqrt(A.xx*A.xx + A.yx*A.yx + A.zx*A.zx),
               sqrt(A.xy*A.xy + A.yy*A.yy + A.zy*A.zy),
               sqrt(A.xz*A.xz + A.yz*A.yz + A.zz*A.zz)};
  double a_0_plus = a_0
    + a.x/(double)gd.x + a.y/(double)gd.y+ a.z/(double)gd.z;
  int nxlim = (int)ceil(as.x*(double)gd.x*a_0_plus - 1);
  int nylim = (int)ceil(as.y*(double)gd.y*a_0_plus - 1);
  int nzlim = (int)ceil(as.z*(double)gd.z*a_0_plus - 1);
  int nxd = 2*nxlim + 1 < gd.x ? 2*nxlim + 1 : gd.x;
  int nyd = 2*nylim + 1 < gd.y ? 2*nylim + 1 : gd.y;
  int nzd = 2*nzlim + 1 < gd.z ? 2*nzlim + 1 : gd.z;
  // compute cell neighbor lists - can be part of preprocessing
  int ndim = nxd*nyd*nzd;
  Triple *n = (Triple *)calloc(ndim, sizeof(Triple));
  int k = 0; // number of occupied elements of n;
  for (int nx = -nxd/2; nx < (nxd + 1)/2; nx++)
    for (int ny = -nyd/2; ny < (nyd + 1)/2; ny++)
      for (int nz = -nzd/2; nz < (nzd + 1)/2; nz++){
        Vector s = {(double)nx/(double)gd.x,
                    (double)ny/(double)gd.y, (double)nz/(double)gd.z};
        Vector r = prod(A, s);
        double normr = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        Vector ATr = prodT(A, r);
        double HATr = fabs(ATr.x)/(double)gd.x
          + fabs(ATr.y)/(double)gd.y + fabs(ATr.z)/(double)gd.z;
        if (normr - HATr/normr < a_0){
          Triple nk = {nx, ny, nz};
          n[k] = nk;
          k++;}
      }
  int kcnt = k;
  // compute interactions
  double beta = ff->beta;
  double beta2 = beta*beta;
  double const_ = beta/sqrt(atan(1.));  // 2 beta/sqrt(pi)
  int m = 0;  // index of grid cell 1
  int i = -1;  // particle number
  while (true){  // loop over i
    if (i >= 0) i = next[i];
    else i = first[m];
    while (i < 0 && m + 1 < gd_prod){m++; i = first[m];}
    //### loop exit ###
    if (i < 0) break;
    int mx = m/(gd.y*gd.z);
    int my = (m - mx*gd.y*gd.z)/gd.z;
    int mz = m - (mx*gd.y + my)*gd.z;
    Vector ri = position[i];
    int k = 0; // index of next neighbor triple for j
    int j = i;
    while (true){  // loop over j
      j = next[j];
      while (j < 0 && k < kcnt){
        int mpnx = (mx + n[k].x + gd.x)%gd.x;
        int mpny = (my + n[k].y + gd.y)%gd.y;
        int mpnz = (mz + n[k].z + gd.z)%gd.z;
        int mpn = (mpnx*gd.y + mpny)*gd.z + mpnz;
        if (mpn > m) j = first[mpn];
        k++;}
      //### loop exit ###
      if (j < 0) break;
      // compute interaction
      Vector rj = position[j];
      Vector r = {rj.x - ri.x, rj.y - ri.y, rj.z - ri.z};
      // convert to nearest image
      Vector s = prod(Ai, r);
      Vector p = {floor(s.x + 0.5), floor(s.y + 0.5), floor(s.z + 0.5)};
      Vector Ap = prod(A, p);
      r.x -= Ap.x; r.y -= Ap.y; r.z -= Ap.z;
      Vector r0 = r;
      double sum = 0.;
      Vector F = {0., 0., 0.};
      for (p.x = - pxlim; p.x <= pxlim; p.x++)
        for (p.y = - pylim; p.y <= pylim; p.y++)
          for (p.z = - pzlim; p.z <= pzlim; p.z++){
            Ap = prod(A, p);
            r.x = r0.x + Ap.x; r.y = r0.y + Ap.y; r.z = r0.z + Ap.z;
            double r2 = r.x*r.x + r.y*r.y + r.z*r.z;
            if (r2 < a_0*a_0){
              // add V(|r|) = erfc(beta|r|)/|r| to sum
              // F on particle i =  V'(|r|)r/|r|
              // = - 1/|r|^2(V(|r|) + 2 beta exp(-beta^2|r|^2)/sqrt(pi)) r
              double r1 = sqrt(r2);
              double Vr1 = erfc(beta*r1)/r1;
              sum += Vr1;
              double coeff = -(Vr1 + const_*exp(-beta2*r2))/r2;
              F.x += coeff*r.x, F.y += coeff*r.y, F.z += coeff*r.z;}}
      energy += charge[i]*charge[j]*sum;
      F.x *= charge[i]*charge[j];
      F.y *= charge[i]*charge[j];
      F.z *= charge[i]*charge[j];
      force[i].x += F.x, force[i].y += F.y,force[i].z += F.z;
      force[j].x -= F.x, force[j].y -= F.y,force[j].z -= F.z;
      // end compute interaction
    } // end loop over j
  } // end loop over i
  free(n);
  free(first);
  free(next);
  return energy;}

static double Bspline(FF *ff, int i, double t);
static void anterpolate(FF *ff, Triple gd, double *q, int N, double *charge,
                        Vector *r){
  Matrix Ai = *(Matrix *)ff->Ai;
  int nu = ff->orderAcc;
  for (int i = 0; i < N; i++){
    Vector ri = r[i];
    Vector s = prod(Ai, ri);
    s.x = s.x - floor(s.x);
    s.y = s.y - floor(s.y);
    s.z = s.z - floor(s.z);
    Vector t = {(double)gd.x*s.x, (double)gd.y*s.y, (double)gd.z*s.z};
    Vector em = {floor(t.x), floor(t.y), floor(t.z)};
    t.x -= em.x, t.y -= em.y, t.z -= em.z;
    Triple m = {(int)em.x, (int)em.y, (int)em.z};
    //**for (int nx = 0; nx < nu; nx++){
      for (int nx = - nu/2; nx < nu/2; nx++){
      //**double Qi = Bspline(ff, nx, t.x);
      double Qi = Bspline(ff, nx + nu/2, t.x);
      //**for (int ny = 0; ny < nu; ny++){
        for (int ny = - nu/2; ny < nu/2; ny++){
        //**double Qj = Bspline(ff, ny, t.y);
        double Qj = Bspline(ff, ny + nu/2, t.y);
        //**for (int nz = 0; nz < nu; nz++){
          for (int nz = - nu/2; nz < nu/2; nz++){
          //**double Qk = Bspline(ff, nz, t.z);
          double Qk = Bspline(ff, nz + nu/2, t.z);
          // add to q_{m+n}
          Triple m_n = {((m.x - nx) % gd.x + gd.x) % gd.x,
                        ((m.y - ny) % gd.y + gd.y) % gd.y,
                        ((m.z - nz) % gd.z + gd.z) % gd.z};
          q[(m_n.x*gd.y + m_n.y)*gd.z + m_n.z] += charge[i]*Qi*Qj*Qk;}}}}}

static void FFTg2g(FF *ff, Triple gd, double *el, double *ql, double *kh){
  // el = DFT of ql
  for (int i = 0; i < gd.x*gd.y*gd.z; i++)
    ff->fftw_in[i] = (fftw_complex)ql[i];
  fftw_execute(ff->forward);
  // ql = kh . el pointwise
  for (int i = 0; i < gd.x*gd.y*gd.z; i++)
        ff->fftw_in[i] *= (fftw_complex)kh[i];
  // el = invDFT of ql
  fftw_execute(ff->backward);
  for (int i = 0; i < gd.x*gd.y*gd.z; i++)
    el[i] = creal(ff->fftw_in[i])/(double)(gd.x*gd.y*gd.z);
}

static double Bspline1(double *Qi1t, FF *ff, int i, double t);
static void interpolate(FF *ff, int N, Vector *E, Vector *r, Triple gd,
 double *el){
  // return grid-level electric field
  Matrix Ai = *(Matrix *)ff->Ai;
  int nu = ff->orderAcc;
  for (int i = 0; i < N; i++){
    double Ex = 0., Ey = 0., Ez = 0.;
    Vector ri = r[i];
    Vector s = prod(Ai, ri);
    Vector t = {(double)gd.x*s.x, (double)gd.y*s.y, (double)gd.z*s.z};
    Vector em = {floor(t.x), floor(t.y), floor(t.z)};
    t.x -= em.x, t.y -= em.y, t.z -= em.z;
    Triple m = {(int)em.x, (int)em.y, (int)em.z};
    //**for (int ni = 0; ni < nu; ni++){
    for (int ni = - nu/2; ni < nu/2; ni++){
      double Qi1;
      //**double Qi = Bspline1(&Qi1, ff, ni, t.x);
      double Qi = Bspline1(&Qi1, ff, ni + nu/2, t.x);
      //**for (int nj = 0; nj < nu; nj++){
      for (int nj = - nu/2; nj < nu/2; nj++){
        double Qj1;
        //**double Qj = Bspline1(&Qj1, ff, nj, t.y);
        double Qj = Bspline1(&Qj1, ff, nj + nu/2, t.y);
        //**for (int nk = 0; nk < nu; nk++){
        for (int nk = - nu/2; nk < nu/2; nk++){
          double Qk1;
          //**double Qk = Bspline1(&Qk1, ff, nk, t.z);
          double Qk = Bspline1(&Qk1, ff, nk + nu/2, t.z);
          Triple m_n = {((m.x - ni) % gd.x + gd.x) % gd.x,
                        ((m.y - nj) % gd.y + gd.y) % gd.y,
                        ((m.z - nk) % gd.z + gd.z) % gd.z};
          double elm_n = el[(m_n.x*gd.y + m_n.y)*gd.z + m_n.z];
          Ex -= Qi1*Qj*Qk*elm_n;
          Ey -= Qi*Qj1*Qk*elm_n;
          Ez -= Qi*Qj*Qk1*elm_n;}}}
    Ex *= (double)gd.x, Ey *= (double)gd.y, Ez *= (double)gd.z;
    Vector E_ = {Ex, Ey, Ez};
    E[i] = prodT(Ai, E_);}}

static double Bspline(FF *ff, int i, double t){
  // assumes 0 <= i < nu and 0 <= t < 1
  int nu = ff->orderAcc;
  if (i >= nu/2){i = nu - 1 - i; t = 1. - t;}
  // piece_i(i + t) = q_i0 + q_i1*t + ... + q_{i,nu-1}*t^{nu-1}
  double *Qi = ff->Q + i*nu;
  double Qit = Qi[nu-1];
  for (int k = nu - 2; k >= 0; k--)
    Qit = Qi[k] + t*Qit;
  return Qit;}

static double Bspline1(double *Qi1t, FF *ff, int i, double t){
  // assumes 0 <= i < nu and 0 <= t < 1
  // piece_i(i + t) = q_i0 + q_i1*t + ... + q_{i,nu-1}*t^{nu-1}
  int nu = ff->orderAcc;
  double sign = 1.;
  if (i >= nu/2){i = nu - 1 - i; t = 1. - t; sign = -1.;}
  double *Qi = ff->Q + i*nu;
  double Qit = Qi[nu-1];
  *Qi1t = 0.;
  for (int k = nu - 2; k >= 0; k--){
    *Qi1t = Qit + t*(*Qi1t);
    Qit = Qi[k] + t*Qit;}
  *Qi1t *= sign;
  return Qit;}

static Vector prod(Matrix M, Vector v){
  Vector Mx = {M.xx*v.x + M.xy*v.y + M.xz*v.z,
               M.yx*v.x + M.yy*v.y + M.yz*v.z,
               M.zx*v.x + M.zy*v.y + M.zz*v.z};
  return Mx;}

static Vector prodT(Matrix M, Vector v){
  Vector MTx = {M.xx*v.x + M.yx*v.y + M.zx*v.z,
                M.xy*v.x + M.yy*v.y + M.zy*v.z,
                M.xz*v.x + M.yz*v.y + M.zz*v.z};
  return MTx;}

