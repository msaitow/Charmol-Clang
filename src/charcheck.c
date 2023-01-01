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
#include "charcheck.h"

// NOTE that changes in charwrite.c file can make changes in this file necessary

/*********************************************************************************************************************************************************************************/
// main functions for checking various parts before writing output
// checking basic settings
void check_basic_settings(int camspec, double *camera, int scenerot, double scenerotangle, int camrot, double *camerarot, int camrq,
                          double *camerarq, int camsky, double *camerasky, int camzoom, double camerazoom, int finish, double *finset)
{
if(camspec==1 && (!isfinite(camera[0]) || !isfinite(camera[1]) || !isfinite(camera[2])))
 {puts("\nWrong number found for cameraposition...\n"); exit(1);}
if(camspec==1 && fabs(camera[0]) < 1E-3 && fabs(camera[1]) < 1E-3 && fabs(camera[2]) < 1E-3)
 {puts("\nCamera position is too close to the origin of coordinate system...\n"); exit(1);}
if(scenerot==1 && !isfinite(scenerotangle))
 {puts("\nWrong number found for scenerotation...\n"); exit(1);}
if(camrot==1 && (!isfinite(camerarot[0]) || !isfinite(camerarot[1]) || !isfinite(camerarot[2])))
 {puts("\nWrong number found for camerarotation...\n"); exit(1);}
if(camrq==1 && (!isfinite(camerarq[0]) || !isfinite(camerarq[1]) || !isfinite(camerarq[2]) || !isfinite(camerarq[3])))
 {puts("\nWrong number found for camerarotquat...\n"); exit(1);}
if(camrq==1 && fabs(camerarq[0]) < 1E-6 && fabs(camerarq[1]) < 1E-6 && fabs(camerarq[2]) < 1E-6 && fabs(camerarq[3]) < 1E-6)
 {puts("\nWrong specification of camerarotquat...\n"); exit(1);}
if(camsky==1 && (!isfinite(camerasky[0]) || !isfinite(camerasky[1]) || !isfinite(camerasky[2])))
 {puts("\nWrong number found for camerasky...\n"); exit(1);}
if(camzoom==1 && !isfinite(camerazoom))
 {puts("\nWrong number found for camerazoom...\n"); exit(1);}
if(finish==1 && (!isfinite(finset[0]) || !isfinite(finset[1]) || !isfinite(finset[2])))
 {puts("\nWrong number found for finish...\n"); exit(1);}
if(!isfinite(finset[3]))
 {puts("\nWrong number found for ambient...\n"); exit(1);}
if(!isfinite(finset[4]))
 {puts("\nWrong number found for shininess...\n"); exit(1);}
}

// checking xyz
void check_xyz(ATOMS *xyz, int natoms, int nmolecules, int hideatoms, int atomcolor, int bondcolor, double *atomradiusscale,
               double *lightatomradiusscale, double *bondradius, double **defatomcolor, double **defbondcolor)
{
int i, j, k;
double rad1, rad2, dist, end[3], rgb[3], bondrad;

// checking atoms
for(i=0;i<natoms;i++)
 {
 if(hideatoms==1)
  {
  rad1 = bondradius[xyz[i].ordid];
  if(bondcolor==1) for(k=0;k<3;k++) rgb[k] = xyz[i].rgb[k];
   else for(k=0;k<3;k++) rgb[k] = defbondcolor[xyz[i].ordid][k];
  }
  else
   {
   if(xyz[i].rad<(0.5*angstrom_to_bohr))
    {
    if(xyz[i].rad*xyz[i].scale*atomradiusscale[xyz[i].ordid]*lightatomradiusscale[xyz[i].ordid] < bondradius[xyz[i].ordid]) rad1 = bondradius[xyz[i].ordid];
     else rad1 = xyz[i].rad*xyz[i].scale*atomradiusscale[xyz[i].ordid]*lightatomradiusscale[xyz[i].ordid];
    }
    else
     {
     if(xyz[i].rad*xyz[i].scale*atomradiusscale[xyz[i].ordid] < bondradius[xyz[i].ordid]) rad1 = bondradius[xyz[i].ordid];
      else rad1 = xyz[i].rad*xyz[i].scale*atomradiusscale[xyz[i].ordid];
     }
   if(atomcolor==1) for(k=0;k<3;k++) rgb[k] = xyz[i].rgb[k];
    else for(k=0;k<3;k++) rgb[k] = defatomcolor[xyz[i].ordid][k];
   }
 if(!isfinite(xyz[i].coor[0]) || !isfinite(xyz[i].coor[1]) || !isfinite(xyz[i].coor[2]))
  {printf("\nWrong number found for position of atom %d...\n\n", i+1); exit(1);}
 if(rad1 < 1E-6 || !isfinite(rad1))
  {printf("\nToo small or wrong number found for radius of atom %d (povray or vrml-player could fail)...\n\n", i+1); exit(1);}
 if(!isfinite(rgb[0]) || !isfinite(rgb[1]) || !isfinite(rgb[2]))
  {printf("\nWrong number found for color of atom %d...\n\n", i+1); exit(1);}
 }
// checking bonds
for(i=0;i<natoms;i++)
 {
 if(hideatoms==1)
  {
  rad1 = bondradius[xyz[i].ordid];
  if(bondcolor==1) for(k=0;k<3;k++) rgb[k] = xyz[i].rgb[k];
   else for(k=0;k<3;k++) rgb[k] = defbondcolor[xyz[i].ordid][k];
  }
  else
   {
   if(xyz[i].rad<(0.5*angstrom_to_bohr))
    {
    if(xyz[i].rad*xyz[i].scale*atomradiusscale[xyz[i].ordid]*lightatomradiusscale[xyz[i].ordid] < bondradius[xyz[i].ordid]) rad1 = bondradius[xyz[i].ordid];
     else rad1 = xyz[i].rad*xyz[i].scale*atomradiusscale[xyz[i].ordid]*lightatomradiusscale[xyz[i].ordid];
    }
    else
     {
     if(xyz[i].rad*xyz[i].scale*atomradiusscale[xyz[i].ordid] < bondradius[xyz[i].ordid]) rad1 = bondradius[xyz[i].ordid];
      else rad1 = xyz[i].rad*xyz[i].scale*atomradiusscale[xyz[i].ordid];
     }
   if(atomcolor==1) for(k=0;k<3;k++) rgb[k] = xyz[i].rgb[k];
    else for(k=0;k<3;k++) rgb[k] = defatomcolor[xyz[i].ordid][k];
   }
 for(j=0;xyz[i].cone[j]!=0;j++)
  {
  if(hideatoms==1) for(k=0;k<3;k++) end[k] = ((xyz[xyz[i].cone[j]-1].coor[k] - xyz[i].coor[k])/2 + xyz[i].coor[k]);
   else
    {
    if(xyz[xyz[i].cone[j]-1].rad<(0.5*angstrom_to_bohr))
     {
     if(xyz[xyz[i].cone[j]-1].rad*xyz[xyz[i].cone[j]-1].scale*atomradiusscale[xyz[xyz[i].cone[j]-1].ordid]*lightatomradiusscale[xyz[xyz[i].cone[j]-1].ordid] < bondradius[xyz[xyz[i].cone[j]-1].ordid])
      rad2 = bondradius[xyz[xyz[i].cone[j]-1].ordid];
      else rad2 = xyz[xyz[i].cone[j]-1].rad*xyz[xyz[i].cone[j]-1].scale*atomradiusscale[xyz[xyz[i].cone[j]-1].ordid]*lightatomradiusscale[xyz[xyz[i].cone[j]-1].ordid];
     }
     else
      {
      if(xyz[xyz[i].cone[j]-1].rad*xyz[xyz[i].cone[j]-1].scale*atomradiusscale[xyz[xyz[i].cone[j]-1].ordid] < bondradius[xyz[xyz[i].cone[j]-1].ordid]) rad2 = bondradius[xyz[xyz[i].cone[j]-1].ordid];
       else rad2 = xyz[xyz[i].cone[j]-1].rad*xyz[xyz[i].cone[j]-1].scale*atomradiusscale[xyz[xyz[i].cone[j]-1].ordid];
      }
    dist = 0.0; for(k=0;k<3;k++) dist += (xyz[xyz[i].cone[j]-1].coor[k] - xyz[i].coor[k])*(xyz[xyz[i].cone[j]-1].coor[k] - xyz[i].coor[k]); dist = sqrt(dist);
    for(k=0;k<3;k++) end[k] = (((dist - rad1 - rad2)/2 + rad1)/dist)*(xyz[xyz[i].cone[j]-1].coor[k] - xyz[i].coor[k]) + xyz[i].coor[k];
    }
  if(bondcolor==1)
   {
   if(atomcolor==1) for(k=0;k<3;k++) rgb[k] = xyz[i].rgb[k];
    else for(k=0;k<3;k++) rgb[k] = defatomcolor[xyz[i].ordid][k];
   }
   else for(k=0;k<3;k++) rgb[k] = defbondcolor[xyz[i].ordid][k];
  if(xyz[i].molid==xyz[xyz[i].cone[j]-1].molid) bondrad = bondradius[xyz[i].ordid];
   else bondrad = bondradius[nmolecules];
  dist = 0.0; for(k=0;k<3;k++) dist += (end[k] - xyz[i].coor[k])*(end[k] - xyz[i].coor[k]); dist = sqrt(dist);
  if(!isfinite(end[0]) || !isfinite(end[1]) || !isfinite(end[2]))
   {printf("\nWrong number found for bond between atoms %d and %d...\n\n", i+1, xyz[i].cone[j]); exit(1);}
  if(dist < 1E-6 || !isfinite(dist))
   {printf("\nToo small or wrong number found for length of bond between atoms %d and %d (povray or vrml-player could fail)...\n\n", i+1, xyz[i].cone[j]); exit(1);}
  if(bondrad < 1E-6 || !isfinite(bondrad))
   {printf("\nToo small or wrong number found for radius of bond between atoms %d and %d (povray or vrml-player could fail)...\n\n", i+1, xyz[i].cone[j]); exit(1);}
  if(!isfinite(rgb[0]) || !isfinite(rgb[1]) || !isfinite(rgb[2]))
   {printf("\nWrong number found for color of bond between atoms %d and %d (povray or vrml-player could fail)...\n\n", i+1, xyz[i].cone[j]); exit(1);}
  }
 }
}

// checking gauges
void check_gauges(ATOMS *xyz, int natoms, GAUGE *gauge, int ngauges)
{
int i, j;
double dist;

for(i=0;i<ngauges;i++)
 {
 dist = 0.0; for(j=0;j<3;j++) dist += (xyz[gauge[i].a2-1].coor[j] - xyz[gauge[i].a1-1].coor[j])*(xyz[gauge[i].a2-1].coor[j] - xyz[gauge[i].a1-1].coor[j]); dist = sqrt(dist);
 if(dist < 1E-6 || !isfinite(dist))
  {printf("\nToo small or wrong number found for length of gauge between atoms %d and %d (povray or vrml-player could fail)...\n\n", gauge[i].a1, gauge[i].a2); exit(1);}
 if(gauge[i].rad < 1E-6 || !isfinite(gauge[i].rad))
  {printf("\nToo small or wrong number found for radius of gauge between atoms %d and %d (povray or vrml-player could fail)...\n\n", gauge[i].a1, gauge[i].a2); exit(1);}
 if(!isfinite(gauge[i].rgb[0]) || !isfinite(gauge[i].rgb[1]) || !isfinite(gauge[i].rgb[2]))
  {printf("\nWrong number found for color of gauge between atoms %d and %d (povray or vrml-player could fail)...\n\n", gauge[i].a1, gauge[i].a2); exit(1);}
 }
}

// checking single arrow
void check_single_arrow(TYCKA tyc, KONUS kon, int arroword, int vib)
{
int i;
double tycdist, kondist;

tycdist = 0.0; for(i=0;i<3;i++) tycdist += (tyc.end[i] - tyc.start[i])*(tyc.end[i] - tyc.start[i]); tycdist = sqrt(tycdist);
kondist = 0.0; for(i=0;i<3;i++) kondist += (kon.end[i] - tyc.end[i])*(kon.end[i] - tyc.end[i]); kondist = sqrt(kondist);

if(tyc.rad < 1E-6 || !isfinite(tyc.rad))
 {
 if(vib==0) {printf("\nToo small or wrong number found for radius of arrow %d (povray or vrml-player could fail)...\n\n", arroword); exit(1);}
  else {printf("\nToo small or wrong number found for radius of arrow for atom %d of vibration %d (povray or vrml-player could fail)...\n\n", arroword, vib); exit(1);}
 }
if(kon.rad < 1E-6 || !isfinite(kon.rad))
 {
 if(vib==0) {printf("\nToo small or wrong number found for radius of arrow %d (povray or vrml-player could fail)...\n\n", arroword); exit(1);}
  else {printf("\nToo small or wrong number found for radius of arrow for atom %d of vibration %d (povray or vrml-player could fail)...\n\n", arroword, vib); exit(1);}
 }
if(tycdist < 1E-6 || !isfinite(tycdist))
 {
 if(vib==0) {printf("\nToo small or wrong number found for length of arrow %d (povray or vrml-player could fail)...\n\n", arroword); exit(1);}
  else {printf("\nToo small or wrong number found for length of arrow for atom %d of vibration %d (povray or vrml-player could fail)...\n\n", arroword, vib); exit(1);}
 }
if(kondist < 1E-6 || !isfinite(kondist))
 {
 if(vib==0) {printf("\nToo small or wrong number found for length of arrow %d (povray or vrml-player could fail)...\n\n", arroword); exit(1);}
  else {printf("\nToo small or wrong number found for length of arrow for atom %d of vibration %d (povray or vrml-player could fail)...\n\n", arroword, vib); exit(1);}
 }
if(!isfinite(tyc.rgb[0]) || !isfinite(tyc.rgb[1]) || !isfinite(tyc.rgb[2]))
 {
 if(vib==0) {printf("\nWrong number found for color of arrow %d (povray or vrml-player could fail)...\n\n", arroword); exit(1);}
  else {printf("\nWrong number found for color of arrow for atom %d of vibration %d (povray or vrml-player could fail)...\n\n", arroword, vib); exit(1);}
 }
}

// checking mesh
void check_mesh(double *rgb, double transparency, double isoval, int surf)
{
if(transparency < 0.0 || !isfinite(transparency))
 {
 if(surf==0) {puts("\nWrong number found for transparency of orbital (povray or vrml-player could fail)...\n"); exit(1);}
  else {puts("\nWrong number found for transparency of surface (povray or vrml-player could fail)...\n"); exit(1);}
 }
if(!isfinite(isoval))
 {
 if(surf==0) {puts("\nWrong number found for isovalue of orbital (povray or vrml-player could fail)...\n"); exit(1);}
  else {puts("\nWrong number found for isovalue of surface (povray or vrml-player could fail)...\n"); exit(1);}
 }
if(!isfinite(rgb[0]) || !isfinite(rgb[1]) || !isfinite(rgb[2]))
 {
 if(surf==0) {puts("\nWrong number found for color of orbital (povray or vrml-player could fail)...\n"); exit(1);}
  else {puts("\nWrong number found for color of surface (povray or vrml-player could fail)...\n"); exit(1);}
 }
}

