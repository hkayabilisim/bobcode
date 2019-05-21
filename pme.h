#include <stdbool.h>
#include <time.h>
#include <fftw3.h>

struct O {
  char *test;
  clock_t time;
}; // struct for optional output
extern struct O o;

// it is recommended that these values be set with the FF_set methods
typedef struct FF {
  int N;
  double A[3][3], Ai[3][3];  // A inverse transpose
  double cutoff;  // abs cutorr
  int orderAcc;
  int topGridDim[3];
  double tolDir;
  double tolRec;
  double kmax;
  fftw_complex *fftw_in;
  fftw_plan forward, backward;
  double beta;
  double *Q;  // B-spline coefficients
  // first ord/2 pieces of B-splines Q_ord(t)
  // piece_i(t) = q_i0 + q_i1*(t - i) + ... + q_{i,ord-1}*(t - i)^{ord-1}
  int kLim[3]; // range of wavenumbers
  double *cL[3]; // c_x^2, c_y^2, c_z^2:
  // coeffs for interpolating reciprocal sum
  double *khat;  // grid2grid stencil; has dimensions topGridDim
  double coeff1, coeff2;  // coeffs of const part of energy
} FF;
FF *FF_new(void);
void FF_set_cutoff(FF *ff, double cutoff);
double FF_get_cutoff(FF *ff);
void FF_set_orderAcc(FF *ff, int orderAcc);
void FF_set_topGridDim(FF *ff, int topGridDim[3]);
void FF_set_tolDir(FF *ff, double tolDir);
void FF_set_tolRec(FF *ff, double tolRec);
void FF_build(FF *ff, int N, double edges[3][3]);
double FF_get_cutoff(FF *ff);
int FF_get_orderAcc(FF *ff);
void FF_get_topGridDim(FF *ff, int topGridDim[3]);
double FF_get_tolDir(FF *ff);
double FF_get_tolRec(FF *ff);
double FF_get_errEst(FF *ff, int N, double *charge);
void FF_rebuild(FF *ff, double edges[3][3]);
double FF_energy(FF *ff, int N, double (*force)[3], double (*position)[3],
                 double *charge, double *weight);
  // if weights == NULL, unit weights are assumed; otherwise
  // weights should point to an array of length FF_get_maxLevel(ff) + 1
void FF_delete(FF *ff);
