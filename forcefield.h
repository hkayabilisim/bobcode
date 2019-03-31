#include <stdbool.h>
#include <time.h>

struct O {
  char *test;
  double *khatLp1;
  double **kappa;
  double e_const;
  double *e; // e[0] = e^0, e[l] = e^{l:}
  bool recX; // e^{L+1} excluded
  clock_t time;
	int level; // 0, 1 (1 to L-1), 2 (L and L+1)
}; // struct for optional output
extern struct O o;

// it is recommended that these values be set with the FF_set methods
typedef struct FF {
  int N;
  double A[3][3], Ai[3][3];  // A inverse transpose
  double relCutoff;
  int orderAcc;
  int maxLevel;
  int topGridDim[3];
  double tolDir;
  double tolRec;
  double kmax;
  double beta;
  double *tau;  // softener coefficients
  double *Q;  // B-spline coefficients
  // first ord/2 pieces of B-splines Q_ord(t)
  // piece_i(t) = q_i0 + q_i1*(t - i) + ... + q_{i,ord-1}*(t - i)^{ord-1}
  double *J;  // 2-scale stencil, indexed from -nu/2 thru nu/2
  double *(*omegap)[3];  // quasi-interpolation coefficients
  int nLim; // ceiling(relCutoff - 1)
  double *aCut; // abs cutoffs
  int kLim[3]; // range of wavenumbers
  double *cL[3]; // c_x^2, c_y^2, c_z^2:
  // coeffs for interpolating reciprocal sum
  double **khat;  // grid2grid stencils
  // khat[L] has dimensions topGridDim
  //khat[l], l < L, has dimensions that are the lesser of  2*nLim + 1
  //  and those of the grid at level l 
  double coeff1, coeff2;  // coeffs of const part of energy
} FF;
FF *FF_new(void);
void FF_set_relCutoff(FF *ff, double rel_cutoff);
void FF_set_orderAcc(FF *ff, int orderAcc);
void FF_set_maxLevel(FF *ff, int maxLevel);
void FF_set_topGridDim(FF *ff, int topGridDim[3]);
void FF_set_tolDir(FF *ff, double tolDir);
void FF_set_tolRec(FF *ff, double tolRec);
void FF_build(FF *ff, int N, double edges[3][3]);
double FF_get_relCutoff(FF *ff);
int FF_get_orderAcc(FF *ff);
int FF_get_maxLevel(FF *ff);
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
