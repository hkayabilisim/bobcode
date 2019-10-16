#include "forcefield.h"
void grid2grid(FF *ff,int l,Triple gd, double *el, double *ql,
                      Triple sd, double *kh){
  double *qlnew = padding_z(ff,l,ql,gd,sd);
  msm4g_tic();  
  int gdznew = gd.z + sd.z ;
  for (int mx = 0; mx < gd.x; mx++) {
    int mxoff = mx * gd.y * gd.z;
    for (int my = 0; my < gd.y; my++) {
      int myoff = mxoff + my * gd.z;
      for (int mz = 0; mz < gd.z; mz++){
        int mzoff = myoff + mz;
        double elsum = 0.0;
        int kxoff = ((mx + sd.x/2 + gd.x)%gd.x) * gd.y * gdznew ;
        for (int nxoff = 0 ; nxoff < sd.x*sd.y*sd.z ; nxoff += sd.y * sd.z){
          int kyoff = kxoff + ((my + sd.y/2 + gd.y)%gd.y) * gdznew ;
          for (int nyoff = nxoff; nyoff < nxoff + sd.y*sd.z ; nyoff += sd.z) {
            int kzoff = kyoff + mz + 2*(sd.z/2);
            for (int nzoff = nyoff; nzoff < nyoff + sd.z; nzoff++){
              elsum += kh[nzoff]*qlnew[kzoff];
              kzoff--;
            }
            kyoff -= gdznew;
            kyoff += (kyoff < kxoff)*gd.y*gdznew;
          }
          kxoff -= gd.y * gdznew;
          kxoff += (kxoff < 0)*gd.x*gd.y*gdznew;
        }
        el[mzoff] += elsum ;
      }
    }
  }
  ff->time_grid2grid[l] = msm4g_toc();
  free(qlnew);
}
