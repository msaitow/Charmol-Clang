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
#include "chararrow.h"
#include "charread.h"
#include "charcheck.h"

/*********************************************************************************************************************************************************************************/
// functions for arrows
// process arrows information and make them
void arrow_make_arrows(TYCKA *tyc, KONUS *kon, int narrows, ATOMS *xyz, int natoms, int nmolecules, double convertxyzval,
                       int centerofmass, double *centerofmasscoor, double *arrowradius, double *coneradiusscale, char **input)
{
int i, j, k, l, c, atomtyc, enlarge;
double dist;

for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "arrow"))!=NULL) break; i++;
for(k=0;k<narrows;k++)
 {
 if((c = sscanf(input[k+i], "%s%lf%lf%lf%lf%lf%lf%lf",
  tyc[k].symb, &tyc[k].scale, &kon[k].end[0], &kon[k].end[1], &kon[k].end[2], &tyc[k].rgb[0], &tyc[k].rgb[1], &tyc[k].rgb[2]))!=8)
   {puts("\nProblems when reading arrows data...\n"); exit(1);}
 for(c=0;c<3;c++) kon[k].end[c] *= convertxyzval;
 if((strstr(tyc[k].symb, "orig"))!=NULL)
  {
  for(c=0;c<3;c++) tyc[k].start[c] = 0.0;
  if(centerofmass==1) {for(c=0;c<3;c++) tyc[k].start[c] -= centerofmasscoor[c];}
  if(nmolecules==1) tyc[k].ordid = 0; else tyc[k].ordid = nmolecules;
  }
  else
   {
   if((strstr(tyc[k].symb, "coor"))!=NULL)
    {
    if((c = sscanf(input[k+i], "%*s%*f%*f%*f%*f%*f%*f%*f%lf%lf%lf",
     &tyc[k].start[0], &tyc[k].start[1], &tyc[k].start[2]))!=3) {puts("\nProblems when reading arrows data...\n"); exit(1);}
    for(c=0;c<3;c++) tyc[k].start[c] *= convertxyzval;
    if(centerofmass==1) {for(c=0;c<3;c++) tyc[k].start[c] -= centerofmasscoor[c];}
    if(nmolecules==1) tyc[k].ordid = 0; else tyc[k].ordid = nmolecules;
    }
    else
     {
     if((c = sscanf(tyc[k].symb, "%d", &atomtyc))!=1) {puts("\nProblems when reading arrows data...\n"); exit(1);}
     if(atomtyc<1 || atomtyc>natoms) {printf("\nAtom number %d out of range for arrow...\n\n", atomtyc); exit(1);}
     for(l=0;l<3;l++) tyc[k].start[l] = xyz[atomtyc-1].coor[l]; tyc[k].ordid = xyz[atomtyc-1].ordid;
     }
   }
 }
// making arrows
enlarge = 0;
for(i=0;i<narrows;i++)
 {
 tyc[i].rad = arrowradius[tyc[i].ordid]; kon[i].rad = 2.0*arrowradius[tyc[i].ordid]*coneradiusscale[tyc[i].ordid];
 for(j=0;j<3;j++) kon[i].end[j] *= tyc[i].scale;
 dist = 0.0; for(j=0;j<3;j++) dist += kon[i].end[j]*kon[i].end[j]; dist = sqrt(dist);
 for(j=0;j<3;j++) tyc[i].end[j] = kon[i].end[j]*(1.0 - kon[i].rad/dist);
 for(j=0;j<3;j++) {tyc[i].end[j] += tyc[i].start[j]; kon[i].end[j] += tyc[i].start[j];}
 // checking for hard problems
 check_single_arrow(tyc[i], kon[i], i+1, 0);
 if(dist > 0.0 && dist < kon[i].rad) enlarge = 1;
 }
// checking for soft problems
if(enlarge==1)
 {
 for(i=0;i<narrows;i++)
  {
  dist = 0.0; for(j=0;j<3;j++) dist += (kon[i].end[j] - tyc[i].start[j])*(kon[i].end[j] - tyc[i].start[j]); dist = sqrt(dist);
  if(dist > 0.0 && dist < kon[i].rad) printf("Enlarge scale for arrow %d in order to display it correctly...\n", i+1);
  }
 }
}

