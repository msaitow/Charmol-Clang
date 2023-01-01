// Copyright (C) 2012 by Jakub Chalupsky (chalupsky.jakub@gmail.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see http://www.gnu.org/licenses/.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "charmol.h"
#include "charomp.h"
#include "charcalc.h"

/*********************************************************************************************************************************************************************************/
// functions for calculation of contributions from AOs
// calculate over primitives in a shell
double calc_prim(double r2, int nprim, double *expon, double *contr)
{
int i;
double val;

val = 0.0;
for(i=0;i<nprim;i++) val += contr[i]*exp(-expon[i]*r2);

return val;
}

// sum over primitives in a shell of generally-contracted basis set
double calc_prim_sum_gridexp(int nprim, double *contr, double *gridexp, int ge)
{
int i;
double val;

val = 0.0;
for(i=0;i<nprim;i++) val += contr[i]*gridexp[ge+i];

return val;
}

// function for calculating screening factor for a shell
double calc_screen_factor_shell(int degen, int nprim, double *expon, double *contr, double *contr2, double r2)
{
double primval, val, x;

x = sqrt(r2);
val = 0.0;

if(degen==1) // orbital s
 {
 primval = calc_prim(r2, nprim, expon, contr);
 val = fabs(primval);
 }
if(degen==3) // orbital p
 {
 primval = calc_prim(r2, nprim, expon, contr);
 val = fabs(primval*x);
 }
if(degen==4) // orbital sp
 {
 primval = calc_prim(r2, nprim, expon, contr);
 val = fabs(primval);
 primval = calc_prim(r2, nprim, expon, contr2);
 if(fabs(primval*x)>val) val = fabs(primval*x);
 }
if(degen==5 || degen==6) // orbital d
 {
 primval = calc_prim(r2, nprim, expon, contr);
 val = fabs(primval*x*x/sqrt(3.0));
 }
if(degen==7 || degen==10) // orbital f
 {
 primval = calc_prim(r2, nprim, expon, contr);
 val = fabs(primval*x*x*x/sqrt(15.0));
 }
if(degen==9 || degen==15) // orbital g
 {
 primval = calc_prim(r2, nprim, expon, contr);
 val = fabs(primval*x*x*x*x/sqrt(105.0));
 }
if(degen==11 || degen==21) // orbital h
 {
 primval = calc_prim(r2, nprim, expon, contr);
 val = fabs(primval*x*x*x*x*x/sqrt(945.0));
 }

return val;
}

// function for predeterming screening factors for the shells based on screenval and their modification for general contraction (taking the biggest value)
void calc_screenr2_values(BASIS *gto, int natoms, int nshell, int maxdegen, double screenval, int gencon)
{
int i, j, l;
double val;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(j, val)
#endif
for(i=0;i<nshell;i++)
 {
 j = 1000;
 while((val = calc_screen_factor_shell(gto[i].degen, gto[i].nprim, gto[i].expon, gto[i].contr, gto[i].contr2, (double) j))>screenval && j <= 10000) j *= 2;
 while((val = calc_screen_factor_shell(gto[i].degen, gto[i].nprim, gto[i].expon, gto[i].contr, gto[i].contr2, (double) j))<screenval && j > 0) j--;
 if(j == 0) gto[i].screenr2 = (double) 10000;
  else gto[i].screenr2 = (double) (j+1);
 }
if(gencon==1)
 {
 for(i=1;i<=natoms;i++)
  {
  for(j=1;j<=maxdegen;j++)
   {
   val = 0.0;
   for(l=0;l<nshell;l++) if(gto[l].atom==i && gto[l].degen==j && gto[l].screenr2>val) val = gto[l].screenr2;
   for(l=0;l<nshell;l++) if(gto[l].atom==i && gto[l].degen==j) gto[l].screenr2 = val;
   }
  }
 }
}

// function for calculating the contribution from a shell for molden format
double calc_shell_molden(double r2, double x, double y, double z, int degen, int nprim, double *expon, double *contr, double *contr2,
                         double *motocalc, int mo, double *gridexp, int ge, int gencon)
{
double val, primval, k1, k2, k3, k4, k5, k6, k7;

val = 0.0;

if(gencon==1) primval = calc_prim_sum_gridexp(nprim, contr, gridexp, ge);
 else primval = calc_prim(r2, nprim, expon, contr);

switch(degen)
 {
 case 1:  // orbital s
  val += motocalc[mo]*primval;                                               // s
  break;
 case 3:  // orbital p (order: px, py, pz)
  val += motocalc[mo]*primval*x; mo++;                                       // px
  val += motocalc[mo]*primval*y; mo++;                                       // py
  val += motocalc[mo]*primval*z;                                             // pz
  break;
 case 4:  // orbital sp (order: s, px, py, pz)
  val += motocalc[mo]*primval; mo++;                                         // s
  if(gencon==1) primval = calc_prim_sum_gridexp(nprim, contr2, gridexp, ge);
   else primval = calc_prim(r2, nprim, expon, contr2);
  val += motocalc[mo]*primval*x; mo++;                                       // px
  val += motocalc[mo]*primval*y; mo++;                                       // py
  val += motocalc[mo]*primval*z;                                             // pz
  break;
 case 5:  // orbital 5d (order: d0, d+1, d-1, d+2, d-2)
  k1 = 1.0/sqrt(12.0); k2 = 1.0/2.0;
  val += motocalc[mo]*primval*(2.0*z*z-x*x-y*y)*k1; mo++;                    // d0
  val += motocalc[mo]*primval*x*z; mo++;                                     // d+1
  val += motocalc[mo]*primval*y*z; mo++;                                     // d-1
  val += motocalc[mo]*primval*(x*x-y*y)*k2; mo++;                            // d+2
  val += motocalc[mo]*primval*x*y;                                           // d-2
  break;
 case 6:  // orbital 6d (order: xx, yy, zz, xy, xz, yz)
  k1 = 1.0/sqrt(3.0);
  val += motocalc[mo]*primval*x*x*k1; mo++;                                  // dxx
  val += motocalc[mo]*primval*y*y*k1; mo++;                                  // dyy
  val += motocalc[mo]*primval*z*z*k1; mo++;                                  // dzz
  val += motocalc[mo]*primval*x*y; mo++;                                     // dxy
  val += motocalc[mo]*primval*x*z; mo++;                                     // dxz
  val += motocalc[mo]*primval*y*z;                                           // dyz
  break;
 case 7:  // orbital 7f (order: f0, f+1, f-1, f+2, f-2, f+3, f-3)
  k1 = 1.0/sqrt(60.0); k2 = 1.0/sqrt(40.0); k3 = 1.0/2.0; k4 = 1.0/sqrt(24.0);
  val += motocalc[mo]*primval*(2.0*z*z*z-3.0*(x*x*z+y*y*z))*k1; mo++;        // f0
  val += motocalc[mo]*primval*(4.0*x*z*z-x*x*x-x*y*y)*k2; mo++;              // f+1
  val += motocalc[mo]*primval*(4.0*y*z*z-y*y*y-x*x*y)*k2; mo++;              // f-1
  val += motocalc[mo]*primval*(x*x*z-y*y*z)*k3; mo++;                        // f+2
  val += motocalc[mo]*primval*x*y*z; mo++;                                   // f-2
  val += motocalc[mo]*primval*(x*x*x-3.0*x*y*y)*k4; mo++;                    // f+3
  val += motocalc[mo]*primval*(3.0*x*x*y-y*y*y)*k4;                          // f-3
  break;
 case 10: // orbital 10f (order: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz)
  k1 = 1.0/sqrt(15.0); k2 = 1.0/sqrt(3.0);
  val += motocalc[mo]*primval*x*x*x*k1; mo++;                                // fxxx
  val += motocalc[mo]*primval*y*y*y*k1; mo++;                                // fyyy
  val += motocalc[mo]*primval*z*z*z*k1; mo++;                                // fzzz
  val += motocalc[mo]*primval*x*y*y*k2; mo++;                                // fxyy
  val += motocalc[mo]*primval*x*x*y*k2; mo++;                                // fxxy
  val += motocalc[mo]*primval*x*x*z*k2; mo++;                                // fxxz
  val += motocalc[mo]*primval*x*z*z*k2; mo++;                                // fxzz
  val += motocalc[mo]*primval*y*z*z*k2; mo++;                                // fyzz
  val += motocalc[mo]*primval*y*y*z*k2; mo++;                                // fyyz
  val += motocalc[mo]*primval*x*y*z;                                         // fxyz
  break;
 case 9:  // orbital 9g (order: g0, g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4)
  k1 = 1.0/sqrt(6720.0); k2 = 1.0/sqrt(168.0); k3 = 1.0/sqrt(336.0); k4 = 1.0/sqrt(84.0);
  k5 = 1.0/sqrt(24.0); k6 = 1.0/sqrt(192.0); k7 = 1.0/sqrt(12.0);
  val += motocalc[mo]*primval*(3.0*x*x*x*x+3.0*y*y*y*y+8.0*z*z*z*z+6.0*x*x*y*y-24.0*(x*x*z*z+y*y*z*z))*k1; mo++;   // g0
  val += motocalc[mo]*primval*(4.0*z*z*z*x-3.0*(y*y*x*z+x*x*x*z))*k2; mo++;                                        // g+1
  val += motocalc[mo]*primval*(4.0*z*z*z*y-3.0*(x*x*y*z+y*y*y*z))*k2; mo++;                                        // g-1
  val += motocalc[mo]*primval*(6.0*x*x*z*z-6.0*y*y*z*z-x*x*x*x+y*y*y*y)*k3; mo++;                                  // g+2
  val += motocalc[mo]*primval*(6.0*z*z*x*y-x*x*x*y-y*y*y*x)*k4; mo++;                                              // g-2
  val += motocalc[mo]*primval*(x*x*x*z-3.0*y*y*x*z)*k5; mo++;                                                      // g+3
  val += motocalc[mo]*primval*(3.0*x*x*y*z-y*y*y*z)*k5; mo++;                                                      // g-3
  val += motocalc[mo]*primval*(x*x*x*x+y*y*y*y-6.0*x*x*y*y)*k6; mo++;                                              // g+4
  val += motocalc[mo]*primval*(x*x*x*y-y*y*y*x)*k7;                                                                // g-4
  break;
 case 15: // orbital 15g (order: xxxx, yyyy, zzzz, xxxy, xxxz, yyyx, yyyz, zzzx, zzzy, xxyy, xxzz, yyzz, xxyz, yyxz, zzxy)
  k1 = 1.0/sqrt(105.0); k2 = 1.0/sqrt(15.0); k3 = 1.0/3.0; k4 = 1.0/sqrt(3.0);
  val += motocalc[mo]*primval*x*x*x*x*k1; mo++;                              // gxxxx
  val += motocalc[mo]*primval*y*y*y*y*k1; mo++;                              // gyyyy
  val += motocalc[mo]*primval*z*z*z*z*k1; mo++;                              // gzzzz
  val += motocalc[mo]*primval*x*x*x*y*k2; mo++;                              // gxxxy
  val += motocalc[mo]*primval*x*x*x*z*k2; mo++;                              // gxxxz
  val += motocalc[mo]*primval*y*y*y*x*k2; mo++;                              // gyyyx
  val += motocalc[mo]*primval*y*y*y*z*k2; mo++;                              // gyyyz
  val += motocalc[mo]*primval*z*z*z*x*k2; mo++;                              // gzzzx
  val += motocalc[mo]*primval*z*z*z*y*k2; mo++;                              // gzzzy
  val += motocalc[mo]*primval*x*x*y*y*k3; mo++;                              // gxxyy
  val += motocalc[mo]*primval*x*x*z*z*k3; mo++;                              // gxxzz
  val += motocalc[mo]*primval*y*y*z*z*k3; mo++;                              // gyyzz
  val += motocalc[mo]*primval*x*x*y*z*k4; mo++;                              // gxxyz
  val += motocalc[mo]*primval*y*y*x*z*k4; mo++;                              // gyyxz
  val += motocalc[mo]*primval*z*z*x*y*k4;                                    // gzzxy
  break;
 }

return val;
}

// function for calculating the contribution from a shell for fchk format
double calc_shell_fchk(double r2, double x, double y, double z, int degen, int nprim, double *expon, double *contr, double *contr2,
                       double *motocalc, int mo, double *gridexp, int ge, int gencon)
{
double val, primval, k1, k2, k3, k4, k5, k6, k7;

val = 0.0;

if(gencon==1) primval = calc_prim_sum_gridexp(nprim, contr, gridexp, ge);
 else primval = calc_prim(r2, nprim, expon, contr);

switch(degen)
 {
 case 1:  // orbital s
  val += motocalc[mo]*primval;                                               // s
  break;
 case 3:  // orbital p (order: px, py, pz)
  val += motocalc[mo]*primval*x; mo++;                                       // px
  val += motocalc[mo]*primval*y; mo++;                                       // py
  val += motocalc[mo]*primval*z;                                             // pz
  break;
 case 4:  // orbital sp (order: s, px, py, pz)
  val += motocalc[mo]*primval; mo++;                                         // s
  if(gencon==1) primval = calc_prim_sum_gridexp(nprim, contr2, gridexp, ge);
   else primval = calc_prim(r2, nprim, expon, contr2);
  val += motocalc[mo]*primval*x; mo++;                                       // px
  val += motocalc[mo]*primval*y; mo++;                                       // py
  val += motocalc[mo]*primval*z;                                             // pz
  break;
 case 5:  // orbital 5d (order: d0, d+1, d-1, d+2, d-2)
  k1 = 1.0/sqrt(12.0); k2 = 1.0/2.0;
  val += motocalc[mo]*primval*(2.0*z*z-x*x-y*y)*k1; mo++;                    // d0
  val += motocalc[mo]*primval*x*z; mo++;                                     // d+1
  val += motocalc[mo]*primval*y*z; mo++;                                     // d-1
  val += motocalc[mo]*primval*(x*x-y*y)*k2; mo++;                            // d+2
  val += motocalc[mo]*primval*x*y;                                           // d-2
  break;
 case 6:  // orbital 6d (order: xx, yy, zz, xy, xz, yz)
  k1 = 1.0/sqrt(3.0);
  val += motocalc[mo]*primval*x*x*k1; mo++;                                  // dxx
  val += motocalc[mo]*primval*y*y*k1; mo++;                                  // dyy
  val += motocalc[mo]*primval*z*z*k1; mo++;                                  // dzz
  val += motocalc[mo]*primval*x*y; mo++;                                     // dxy
  val += motocalc[mo]*primval*x*z; mo++;                                     // dxz
  val += motocalc[mo]*primval*y*z;                                           // dyz
  break;
 case 7:  // orbital 7f (order: f0, f+1, f-1, f+2, f-2, f+3, f-3)
  k1 = 1.0/sqrt(60.0); k2 = 1.0/sqrt(40.0); k3 = 1.0/2.0; k4 = 1.0/sqrt(24.0);
  val += motocalc[mo]*primval*(2.0*z*z*z-3.0*(x*x*z+y*y*z))*k1; mo++;        // f0
  val += motocalc[mo]*primval*(4.0*x*z*z-x*x*x-x*y*y)*k2; mo++;              // f+1
  val += motocalc[mo]*primval*(4.0*y*z*z-y*y*y-x*x*y)*k2; mo++;              // f-1
  val += motocalc[mo]*primval*(x*x*z-y*y*z)*k3; mo++;                        // f+2
  val += motocalc[mo]*primval*x*y*z; mo++;                                   // f-2
  val += motocalc[mo]*primval*(x*x*x-3.0*x*y*y)*k4; mo++;                    // f+3
  val += motocalc[mo]*primval*(3.0*x*x*y-y*y*y)*k4;                          // f-3
  break;
 case 10: // orbital 10f (order: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz)
  k1 = 1.0/sqrt(15.0); k2 = 1.0/sqrt(3.0);
  val += motocalc[mo]*primval*x*x*x*k1; mo++;                                // fxxx
  val += motocalc[mo]*primval*y*y*y*k1; mo++;                                // fyyy
  val += motocalc[mo]*primval*z*z*z*k1; mo++;                                // fzzz
  val += motocalc[mo]*primval*x*y*y*k2; mo++;                                // fxyy
  val += motocalc[mo]*primval*x*x*y*k2; mo++;                                // fxxy
  val += motocalc[mo]*primval*x*x*z*k2; mo++;                                // fxxz
  val += motocalc[mo]*primval*x*z*z*k2; mo++;                                // fxzz
  val += motocalc[mo]*primval*y*z*z*k2; mo++;                                // fyzz
  val += motocalc[mo]*primval*y*y*z*k2; mo++;                                // fyyz
  val += motocalc[mo]*primval*x*y*z;                                         // fxyz
  break;
 case 9:  // orbital 9g (order: g0, g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4)
  k1 = 1.0/sqrt(6720.0); k2 = 1.0/sqrt(168.0); k3 = 1.0/sqrt(336.0); k4 = 1.0/sqrt(84.0);
  k5 = 1.0/sqrt(24.0); k6 = 1.0/sqrt(192.0); k7 = 1.0/sqrt(12.0);
  val += motocalc[mo]*primval*(3.0*x*x*x*x+3.0*y*y*y*y+8.0*z*z*z*z+6.0*x*x*y*y-24.0*(x*x*z*z+y*y*z*z))*k1; mo++;   // g0
  val += motocalc[mo]*primval*(4.0*z*z*z*x-3.0*(y*y*x*z+x*x*x*z))*k2; mo++;                                        // g+1
  val += motocalc[mo]*primval*(4.0*z*z*z*y-3.0*(x*x*y*z+y*y*y*z))*k2; mo++;                                        // g-1
  val += motocalc[mo]*primval*(6.0*x*x*z*z-6.0*y*y*z*z-x*x*x*x+y*y*y*y)*k3; mo++;                                  // g+2
  val += motocalc[mo]*primval*(6.0*z*z*x*y-x*x*x*y-y*y*y*x)*k4; mo++;                                              // g-2
  val += motocalc[mo]*primval*(x*x*x*z-3.0*y*y*x*z)*k5; mo++;                                                      // g+3
  val += motocalc[mo]*primval*(3.0*x*x*y*z-y*y*y*z)*k5; mo++;                                                      // g-3
  val += motocalc[mo]*primval*(x*x*x*x+y*y*y*y-6.0*x*x*y*y)*k6; mo++;                                              // g+4
  val += motocalc[mo]*primval*(x*x*x*y-y*y*y*x)*k7;                                                                // g-4
  break;
 case 15: // orbital 15g (order: zzzz, yzzz, yyzz, yyyz, yyyy, xzzz, xyzz, xyyz, xyyy, xxzz, xxyz, xxyy, xxxz, xxxy, xxxx)
  k1 = 1.0/sqrt(105.0); k2 = 1.0/sqrt(15.0); k3 = 1.0/3.0; k4 = 1.0/sqrt(3.0);
  val += motocalc[mo]*primval*z*z*z*z*k1; mo++;                              // gzzzz
  val += motocalc[mo]*primval*z*z*z*y*k2; mo++;                              // gzzzy
  val += motocalc[mo]*primval*y*y*z*z*k3; mo++;                              // gyyzz
  val += motocalc[mo]*primval*y*y*y*z*k2; mo++;                              // gyyyz
  val += motocalc[mo]*primval*y*y*y*y*k1; mo++;                              // gyyyy
  val += motocalc[mo]*primval*z*z*z*x*k2; mo++;                              // gzzzx
  val += motocalc[mo]*primval*z*z*x*y*k4; mo++;                              // gzzxy
  val += motocalc[mo]*primval*y*y*x*z*k4; mo++;                              // gyyxz
  val += motocalc[mo]*primval*y*y*y*x*k2; mo++;                              // gyyyx
  val += motocalc[mo]*primval*x*x*z*z*k3; mo++;                              // gxxzz
  val += motocalc[mo]*primval*x*x*y*z*k4; mo++;                              // gxxyz
  val += motocalc[mo]*primval*x*x*y*y*k3; mo++;                              // gxxyy
  val += motocalc[mo]*primval*x*x*x*z*k2; mo++;                              // gxxxz
  val += motocalc[mo]*primval*x*x*x*y*k2; mo++;                              // gxxxy
  val += motocalc[mo]*primval*x*x*x*x*k1;                                    // gxxxx
  break;
 case 11: // orbital 11h (order: h0, h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5)
  k1 = 1.0/(8.0*sqrt(945.0)); k2 = sqrt(15.0)/(8.0*sqrt(945.0)); k3 = sqrt(105.0)/(4.0*sqrt(945.0));
  k4 = sqrt(70.0)/(16.0*sqrt(945.0)); k5 = 3.0*sqrt(35.0)/(8.0*sqrt(945.0)); k6 = 3.0*sqrt(14.0)/(16.0*sqrt(945.0));
  val += motocalc[mo]*primval*z*(63.0*z*z*z*z - 70.0*z*z*r2 + 15.0*r2*r2)*k1;   mo++;                              // h0
  val += motocalc[mo]*primval*(21.0*z*z*z*z - 14.0*z*z*r2 + r2*r2)*x*k2;        mo++;                              // h+1
  val += motocalc[mo]*primval*(21.0*z*z*z*z - 14.0*z*z*r2 + r2*r2)*y*k2;        mo++;                              // h-1
  val += motocalc[mo]*primval*z*(3.0*z*z - r2)*(x*x - y*y)*k3;                  mo++;                              // h+2
  val += motocalc[mo]*primval*z*(3.0*z*z - r2)*(2.0*x*y)*k3;                    mo++;                              // h-2
  val += motocalc[mo]*primval*(9.0*z*z - r2)*(x*x*x - 3.0*x*y*y)*k4;            mo++;                              // h+3
  val += motocalc[mo]*primval*(9.0*z*z - r2)*(3.0*x*x*y - y*y*y)*k4;            mo++;                              // h-3
  val += motocalc[mo]*primval*z*(x*x*x*x - 6.0*x*x*y*y + y*y*y*y)*k5;           mo++;                              // h+4
  val += motocalc[mo]*primval*z*(4.0*x*x*x*y - 4.0*x*y*y*y)*k5;                 mo++;                              // h-4
  val += motocalc[mo]*primval*(x*x*x*x*x - 10.0*x*x*x*y*y + 5.0*x*y*y*y*y)*k6;  mo++;                              // h+5
  val += motocalc[mo]*primval*(5.0*x*x*x*x*y - 10.0*x*x*y*y*y + y*y*y*y*y)*k6;                                     // h-5
  break;
 case 21: // orbital 21h (order: zzzzz, yzzzz, yyzzz, yyyzz, yyyyz, yyyyy, xzzzz, xyzzz, xyyzz, xyyyz, xyyyy, xxzzz, xxyzz, xxyyz, xxyyy, xxxzz, xxxyz, xxxyy, xxxxz, xxxxy, xxxxx)
  k1 = 1.0/sqrt(945.0); k2 = 1.0/sqrt(105.0); k3 = 1.0/sqrt(45.0); k4 = 1.0/sqrt(15.0); k5 = 1.0/3.0;
  val += motocalc[mo]*primval*z*z*z*z*z*k1; mo++;                            // hzzzzz
  val += motocalc[mo]*primval*y*z*z*z*z*k2; mo++;                            // hyzzzz
  val += motocalc[mo]*primval*y*y*z*z*z*k3; mo++;                            // hyyzzz
  val += motocalc[mo]*primval*y*y*y*z*z*k3; mo++;                            // hyyyzz
  val += motocalc[mo]*primval*y*y*y*y*z*k2; mo++;                            // hyyyyz
  val += motocalc[mo]*primval*y*y*y*y*y*k1; mo++;                            // hyyyyy
  val += motocalc[mo]*primval*x*z*z*z*z*k2; mo++;                            // hxzzzz
  val += motocalc[mo]*primval*x*y*z*z*z*k4; mo++;                            // hxyzzz
  val += motocalc[mo]*primval*x*y*y*z*z*k5; mo++;                            // hxyyzz
  val += motocalc[mo]*primval*x*y*y*y*z*k4; mo++;                            // hxyyyz
  val += motocalc[mo]*primval*x*y*y*y*y*k2; mo++;                            // hxyyyy
  val += motocalc[mo]*primval*x*x*z*z*z*k3; mo++;                            // hxxzzz
  val += motocalc[mo]*primval*x*x*y*z*z*k5; mo++;                            // hxxyzz
  val += motocalc[mo]*primval*x*x*y*y*z*k5; mo++;                            // hxxyyz
  val += motocalc[mo]*primval*x*x*y*y*y*k3; mo++;                            // hxxyyy
  val += motocalc[mo]*primval*x*x*x*z*z*k3; mo++;                            // hxxxzz
  val += motocalc[mo]*primval*x*x*x*y*z*k4; mo++;                            // hxxxyz
  val += motocalc[mo]*primval*x*x*x*y*y*k3; mo++;                            // hxxxyy
  val += motocalc[mo]*primval*x*x*x*x*z*k2; mo++;                            // hxxxxz
  val += motocalc[mo]*primval*x*x*x*x*y*k2; mo++;                            // hxxxxy
  val += motocalc[mo]*primval*x*x*x*x*x*k1;                                  // hxxxxx
  break;
 }

return val;
}

/*********************************************************************************************************************************************************************************/
// main functions for orbital calculations

// with general contraction
// for molden
void calc_grid_gencon_molden(BASIS *gto, int natoms, int nshell, int maxdegen, int totpoints, int screen, double *gridx, double *gridy, double *gridz, double *gridval,
                             double *motocalc, double *gridexp, int gencon)
{
int i, j, k, l, m, c, atom, degen, nprim;
double gx, gy, gz, ax, ay, az, r2, *expon, *contr, *contr2, screenr2;

for(i=0;i<totpoints;i++) gridval[i] = 0.0;
for(i=1;i<=natoms;i++)
 {
 for(j=1;j<=maxdegen;j++)
  {
  c = 0;
  for(k=0;k<nshell;k++) {if(gto[k].atom==i && gto[k].degen==j) break; else c += gto[k].degen;}
  if(k<nshell)
   {
   atom = gto[k].atom;
   degen = gto[k].degen;
   nprim = gto[k].nprim;
   expon = gto[k].expon;
   screenr2 = gto[k].screenr2;
   gx = gto[k].coor[0];
   gy = gto[k].coor[1];
   gz = gto[k].coor[2];
#ifdef _OPENMP
   #pragma omp parallel for default(shared) private(ax, ay, az, r2, m)
#endif
   for(l=0;l<totpoints;l++)
    {
    ax = gridx[l] - gx;
    ay = gridy[l] - gy;
    az = gridz[l] - gz;
    r2 = ax*ax + ay*ay + az*az;
    if(screen==1) {if(r2<screenr2 && r2>1E-10) for(m=0;m<nprim;m++) gridexp[l*nprim+m] = exp(-expon[m]*r2);}
     else {if(r2>1E-10) for(m=0;m<nprim;m++) gridexp[l*nprim+m] = exp(-expon[m]*r2);}
    }
   while(k<nshell)
    {
    if(gto[k].atom==atom && gto[k].degen==degen)
     {
     contr = gto[k].contr;
     contr2 = gto[k].contr2;
#ifdef _OPENMP
     #pragma omp parallel for default(shared) private(ax, ay, az, r2)
#endif
     for(l=0;l<totpoints;l++)
      {
      ax = gridx[l] - gx;
      ay = gridy[l] - gy;
      az = gridz[l] - gz;
      r2 = ax*ax + ay*ay + az*az;
      if(screen==1) {if(r2<screenr2 && r2>1E-10) gridval[l] += calc_shell_molden(r2, ax, ay, az, degen, nprim, expon, contr, contr2, motocalc, c, gridexp, l*nprim, gencon);}
       else {if(r2>1E-10) gridval[l] += calc_shell_molden(r2, ax, ay, az, degen, nprim, expon, contr, contr2, motocalc, c, gridexp, l*nprim, gencon);}
      }
     }
    c += gto[k].degen;
    k++;
    }
   }
  }
 }
}

// for fchk
void calc_grid_gencon_fchk(BASIS *gto, int natoms, int nshell, int maxdegen, int totpoints, int screen, double *gridx, double *gridy, double *gridz, double *gridval,
                           double *motocalc, double *gridexp, int gencon)
{
int i, j, k, l, m, c, atom, degen, nprim;
double gx, gy, gz, ax, ay, az, r2, *expon, *contr, *contr2, screenr2;

for(i=0;i<totpoints;i++) gridval[i] = 0.0;
for(i=1;i<=natoms;i++)
 {
 for(j=1;j<=maxdegen;j++)
  {
  c = 0;
  for(k=0;k<nshell;k++) {if(gto[k].atom==i && gto[k].degen==j) break; else c += gto[k].degen;}
  if(k<nshell)
   {
   atom = gto[k].atom;
   degen = gto[k].degen;
   nprim = gto[k].nprim;
   expon = gto[k].expon;
   screenr2 = gto[k].screenr2;
   gx = gto[k].coor[0];
   gy = gto[k].coor[1];
   gz = gto[k].coor[2];
#ifdef _OPENMP
   #pragma omp parallel for default(shared) private(ax, ay, az, r2, m)
#endif
   for(l=0;l<totpoints;l++)
    {
    ax = gridx[l] - gx;
    ay = gridy[l] - gy;
    az = gridz[l] - gz;
    r2 = ax*ax + ay*ay + az*az;
    if(screen==1) {if(r2<screenr2 && r2>1E-10) for(m=0;m<nprim;m++) gridexp[l*nprim+m] = exp(-expon[m]*r2);}
     else {if(r2>1E-10) for(m=0;m<nprim;m++) gridexp[l*nprim+m] = exp(-expon[m]*r2);}
    }
   while(k<nshell)
    {
    if(gto[k].atom==atom && gto[k].degen==degen)
     {
     contr = gto[k].contr;
     contr2 = gto[k].contr2;
#ifdef _OPENMP
     #pragma omp parallel for default(shared) private(ax, ay, az, r2)
#endif
     for(l=0;l<totpoints;l++)
      {
      ax = gridx[l] - gx;
      ay = gridy[l] - gy;
      az = gridz[l] - gz;
      r2 = ax*ax + ay*ay + az*az;
      if(screen==1) {if(r2<screenr2 && r2>1E-10) gridval[l] += calc_shell_fchk(r2, ax, ay, az, degen, nprim, expon, contr, contr2, motocalc, c, gridexp, l*nprim, gencon);}
       else {if(r2>1E-10) gridval[l] += calc_shell_fchk(r2, ax, ay, az, degen, nprim, expon, contr, contr2, motocalc, c, gridexp, l*nprim, gencon);}
      }
     }
    c += gto[k].degen;
    k++;
    }
   }
  }
 }
}

// without general contraction
// for molden
void calc_grid_segcon_molden(BASIS *gto, int natoms, int nshell, int maxdegen, int totpoints, int screen, double *gridx, double *gridy, double *gridz, double *gridval,
                             double *motocalc, double *gridexp, int gencon)
{
int i, j, c, degen, nprim;
double gx, gy, gz, ax, ay, az, r2, val, *expon, *contr, *contr2, screenr2;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(j, gx, gy, gz, ax, ay, az, r2, val, c, degen, nprim, expon, contr, contr2, screenr2)
#endif
for(i=0;i<totpoints;i++)
 {
 gx = gridx[i]; gy = gridy[i]; gz = gridz[i]; val = 0.0; c = 0;
 for(j=0;j<nshell;j++)
  {
  degen = gto[j].degen;
  nprim = gto[j].nprim;
  expon = gto[j].expon;
  contr = gto[j].contr;
  contr2 = gto[j].contr2;
  screenr2 = gto[j].screenr2;
  ax = gx - gto[j].coor[0];
  ay = gy - gto[j].coor[1];
  az = gz - gto[j].coor[2];
  r2 = ax*ax + ay*ay + az*az;
  if(screen==1) {if(r2<screenr2 && r2>1E-10) val += calc_shell_molden(r2, ax, ay, az, degen, nprim, expon, contr, contr2, motocalc, c, gridexp, 0, gencon);}
   else {if(r2>1E-10) val += calc_shell_molden(r2, ax, ay, az, degen, nprim, expon, contr, contr2, motocalc, c, gridexp, 0, gencon);}
  c += degen;
  }
 gridval[i] = val;
 }
}

// for fchk
void calc_grid_segcon_fchk(BASIS *gto, int natoms, int nshell, int maxdegen, int totpoints, int screen, double *gridx, double *gridy, double *gridz, double *gridval,
                           double *motocalc, double *gridexp, int gencon)
{
int i, j, c, degen, nprim;
double gx, gy, gz, ax, ay, az, r2, val, *expon, *contr, *contr2, screenr2;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(j, gx, gy, gz, ax, ay, az, r2, val, c, degen, nprim, expon, contr, contr2, screenr2)
#endif
for(i=0;i<totpoints;i++)
 {
 gx = gridx[i]; gy = gridy[i]; gz = gridz[i]; val = 0.0; c = 0;
 for(j=0;j<nshell;j++)
  {
  degen = gto[j].degen;
  nprim = gto[j].nprim;
  expon = gto[j].expon;
  contr = gto[j].contr;
  contr2 = gto[j].contr2;
  screenr2 = gto[j].screenr2;
  ax = gx - gto[j].coor[0];
  ay = gy - gto[j].coor[1];
  az = gz - gto[j].coor[2];
  r2 = ax*ax + ay*ay + az*az;
  if(screen==1) {if(r2<screenr2 && r2>1E-10) val += calc_shell_fchk(r2, ax, ay, az, degen, nprim, expon, contr, contr2, motocalc, c, gridexp, 0, gencon);}
   else {if(r2>1E-10) val += calc_shell_fchk(r2, ax, ay, az, degen, nprim, expon, contr, contr2, motocalc, c, gridexp, 0, gencon);}
  c += degen;
  }
 gridval[i] = val;
 }
}

/*********************************************************************************************************************************************************************************/
// other functions connected to orbitals and basis set
// function for checking whether it is generally-contracted basis set
int calc_check_gencon(int natoms, int nshell, int maxdegen, BASIS *gto)
{
int i, j, k, l, c, val;

val = 1;
for(i=1;i<=natoms && val==1;i++)
 {
 for(j=1;j<=maxdegen && val==1;j++)
  {
  for(l=0;l<nshell;l++) if(gto[l].atom==i && gto[l].degen==j) break;
  if(l<nshell)
   {
   k = l+1;
   while(k<nshell)
    {
    if(gto[k].atom==i && gto[k].degen==j)
     {
     if(gto[k].nprim==gto[l].nprim)
      {
      for(c=0;c<gto[k].nprim;c++) if(fabs(gto[k].expon[c]-gto[l].expon[c])>1E-10) val = 0;
      }
      else val = 0;
     }
    k++;
    }
   }
  }
 }

return val;
}

