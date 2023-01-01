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
#include "charatoms.h"
#include "charread.h"

/*********************************************************************************************************************************************************************************/
// functions for xyz units
// determining units of xyz geometry from various inputs
double atoms_check_xyz_units(int extfileformat, char **input, char **extinput)
{
int i, j, c;
double val = 1.0;

switch(extfileformat)
 {
 case 0: // charmol input (default angstroms, units can be specified)
  for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "xyz"))!=NULL) break;
  if((strstr(input[i], "au"))!=NULL) val = 1.0;
   else val = angstrom_to_bohr;
  break;
 case 1: // molden (ATOMS - units specified, simple - angstroms, FR-COORD - atomic units)
  for(i=0;extinput[i]!=NULL;i++) if((strstr(extinput[i], "[ATOMS]"))!=NULL || (strstr(extinput[i], "[Atoms]"))!=NULL || (strstr(extinput[i], "[atoms]"))!=NULL) break;
  if(extinput[i]!=NULL)
   {
   if((strstr(extinput[i], "ANGS"))!=NULL || (strstr(extinput[i], "Angs"))!=NULL || (strstr(extinput[i], "angs"))!=NULL) val = angstrom_to_bohr;
    else val = 1.0;
   }
   else
    {
    if((c = sscanf(extinput[0], "%d", &j))==1) val = angstrom_to_bohr;
     else
      {
      for(i=0;extinput[i]!=NULL;i++) if((strstr(extinput[i], "[FR-COORD]"))!=NULL || (strstr(extinput[i], "[Fr-coord]"))!=NULL || (strstr(extinput[i], "[fr-coord]"))!=NULL) break;
      if(extinput[i]!=NULL) val = 1.0;
      }
    }
  break;
 case 2: // fchk (xyz in atomic units)
  val = 1.0;
  break;
 case 3: // cube (xyz in atomic units)
  val = 1.0;
  break;
 }

return val;
}

// rescale xyz geometry
void atoms_rescale_xyz(ATOMS *xyz, int natoms, double scaleval)
{
int i, j;

for(i=0;i<natoms;i++)
 {
 for(j=0;j<3;j++) xyz[i].coor[j] *= scaleval;
 }
}

/*********************************************************************************************************************************************************************************/
// functions for settings of atoms
// fill basic settings
void atoms_fill_basic_settings(ATOMS *xyz, int natoms, char **atoms)
{
int i, j, c, numtest;
char test[10];

for(i=0;i<natoms;i++)
 {
 if(xyz[i].symb[0]!='\0')
  {
  if((strstr(xyz[i].symb, "Center:"))!=NULL)
   {
   j = 0; while(xyz[i].symb[j]!=':') j++;
   for(c=0,j++;xyz[i].symb[j]!='\0';c++,j++) xyz[i].symb[c] = xyz[i].symb[j];
   xyz[i].symb[c] = xyz[i].symb[j];
   }
  if(xyz[i].symb[0] >= 'a' && xyz[i].symb[0] <= 'z') xyz[i].symb[0] = xyz[i].symb[0] - 'a' + 'A';
  if(xyz[i].symb[1] >= 'A' && xyz[i].symb[1] <= 'Z') xyz[i].symb[1] = xyz[i].symb[1] - 'A' + 'a';
  if(xyz[i].symb[1] >= '0' && xyz[i].symb[1] <= '9') xyz[i].symb[1] = '\0';
  if(xyz[i].symb[2] >= '0' && xyz[i].symb[2] <= '9') xyz[i].symb[2] = '\0';
  strcpy(test, xyz[i].symb);
  strcat(test, " ");
  for(j=0;atoms[j]!=NULL;j++)
   {
   if((strstr(atoms[j], test))!=NULL)
    {
    if((c = sscanf(atoms[j], "%*s%lf%lf%lf%lf%lf%d", &xyz[i].rad, &xyz[i].rgb[0], &xyz[i].rgb[1], &xyz[i].rgb[2], &xyz[i].mass, &xyz[i].atomnum))!=6)
     {puts("\nProblems when reading atoms info...\n"); exit(1);}
    break;
    }
   }
  if(atoms[j]==NULL)
   {printf("\nSpecification for atom %s not set! Check your atomsinfofile or add entry to charatomsinfo.c and recompile...\n\n", xyz[i].symb); exit(1);}
  }
  else
   {
   for(j=0;atoms[j]!=NULL;j++)
    {
    if((c = sscanf(atoms[j], "%*s%*f%*f%*f%*f%*f%d", &numtest))!=1) {puts("\nProblems when reading atoms info...\n"); exit(1);}
    if(numtest==xyz[i].atomnum)
     {
     if((c = sscanf(atoms[j], "%s%lf%lf%lf%lf%lf", xyz[i].symb, &xyz[i].rad, &xyz[i].rgb[0], &xyz[i].rgb[1], &xyz[i].rgb[2], &xyz[i].mass))!=6)
      {puts("\nProblems when reading atoms info...\n"); exit(1);}
     break;
     }
    }
   if(atoms[j]==NULL)
    {printf("\nSpecification for atomic number %d not set! Check your atomsinfofile or add entry to charatomsinfo.c and recompile...\n\n", xyz[i].atomnum); exit(1);}
   }
 xyz[i].scale = 1.0;
 xyz[i].rad *= angstrom_to_bohr;
 }
}

// user-specified changes of atomic scales
void atoms_change_atomscales(ATOMS *xyz, int natoms, char **input)
{
int i, j, k, c, natomscale, atomnum;
char atomsymb[10];
double atomscale;

for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "atomscale"))!=NULL) break; i++;
for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endatomscale"))!=NULL) break;
if(input[i-1]==NULL) natomscale = 0;
 else natomscale = j - i;
if(natomscale!=0)
 {
 for(k=0;k+i<j;k++)
  {
  if((c = sscanf(input[k+i], "%s%lf", atomsymb, &atomscale))!=2) {puts("\nProblems when reading atomscale data...\n"); exit(1);}
  if(atomsymb[0] >= '0' && atomsymb[0] <= '9')
   {
   if((c = sscanf(atomsymb, "%d", &atomnum))!=1) {puts("\nProblems when reading atomscale data...\n"); exit(1);}
   if(atomnum<1 || atomnum>natoms) {printf("\nAtom number %d out of range for atomscale...\n\n", atomnum); exit(1);}
   xyz[atomnum-1].scale = atomscale;
   }
   else
    {
    if(atomsymb[0] >= 'a' && atomsymb[0] <= 'z') atomsymb[0] = atomsymb[0] - 'a' + 'A';
    if(atomsymb[1] >= 'A' && atomsymb[1] <= 'Z') atomsymb[1] = atomsymb[1] - 'A' + 'a';
    if(atomsymb[1] >= '0' && atomsymb[1] <= '9') atomsymb[1] = '\0';
    if(atomsymb[2] >= '0' && atomsymb[2] <= '9') atomsymb[2] = '\0';
    for(c=0;c<natoms;c++)
     {
     if((strcmp(atomsymb, xyz[c].symb))==0) xyz[c].scale = atomscale;
     }
    }
  }
 }
}

// user-specified changes of atomic colors
void atoms_change_atomcolors(ATOMS *xyz, int natoms, char **input)
{
int i, j, k, c, sl, natomcolor, atomnum;
char atomsymb[10];
double atomcolorspec[3];

for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "atomcolor"))!=NULL) break; i++;
for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endatomcolor"))!=NULL) break;
if(input[i-1]==NULL) natomcolor = 0;
 else natomcolor = j - i;
if(natomcolor!=0)
 {
 for(k=0;k+i<j;k++)
  {
  if((c = sscanf(input[k+i], "%s%lf%lf%lf", atomsymb, &atomcolorspec[0], &atomcolorspec[1], &atomcolorspec[2]))!=4) {puts("\nProblems when reading atomcolor data...\n"); exit(1);}
  if(atomsymb[0] >= '0' && atomsymb[0] <= '9')
   {
   if((c = sscanf(atomsymb, "%d", &atomnum))!=1) {puts("\nProblems when reading atomcolor data...\n"); exit(1);}
   if(atomnum<1 || atomnum>natoms) {printf("\nAtom number %d out of range for atomcolor...\n\n", atomnum); exit(1);}
   for(c=0;c<3;c++) xyz[atomnum-1].rgb[c] = atomcolorspec[c];
   }
   else
    {
    if(atomsymb[0] >= 'a' && atomsymb[0] <= 'z') atomsymb[0] = atomsymb[0] - 'a' + 'A';
    if(atomsymb[1] >= 'A' && atomsymb[1] <= 'Z') atomsymb[1] = atomsymb[1] - 'A' + 'a';
    if(atomsymb[1] >= '0' && atomsymb[1] <= '9') atomsymb[1] = '\0';
    if(atomsymb[2] >= '0' && atomsymb[2] <= '9') atomsymb[2] = '\0';
    for(c=0;c<natoms;c++)
     {
     if((strcmp(atomsymb, xyz[c].symb))==0) {for(sl=0;sl<3;sl++) xyz[c].rgb[sl] = atomcolorspec[sl];}
     }
    }
  }
 }
}

/*********************************************************************************************************************************************************************************/
// functions for connectivity matrix
// making connectivity matrix and applying modifications specified by user
void atoms_make_connectivitymatrix(ATOMS *xyz, int natoms, char **input)
{
int i, j, k, sl, c, atom1, atom2;
double dist;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(sl, j, dist, k)
#endif
for(i=0;i<natoms;i++)
 {
 sl = 0;
 if((xyz[i].cone = (int *) malloc((sl+1)*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 for(j=0;j<natoms;j++)
  {
  if(xyz[i].molid == xyz[j].molid)
   {
   dist = 0.0; for(k=0;k<3;k++) dist += (xyz[i].coor[k] - xyz[j].coor[k])*(xyz[i].coor[k] - xyz[j].coor[k]); dist = sqrt(dist);
   if(i!=j && dist<1.2*(xyz[i].rad + xyz[j].rad))
    {
    xyz[i].cone[sl] = j+1;
    sl++;
    if((xyz[i].cone = (int *) realloc(xyz[i].cone, (sl+1)*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
    }
   }
  }
 xyz[i].cone[sl] = 0;
 }

for(i=0;input[i]!=NULL;i++)
 {
 if((read_firststrcmp(input[i], "addbond"))!=NULL)
  {
  if((c = sscanf(input[i], "%*s%d%d", &atom1, &atom2))!=2) {puts("\nProblems when reading addbond data...\n"); exit(1);}
  if(atom1!=atom2)
   {
   if(atom1<1 || atom1>natoms) {printf("\nAtom number %d out of range for addbond...\n\n", atom1); exit(1);}
   if(atom2<1 || atom2>natoms) {printf("\nAtom number %d out of range for addbond...\n\n", atom2); exit(1);}
   j=atom1-1; k=0; while(xyz[j].cone[k]!=0) k++;
   if((xyz[j].cone = (int *) realloc(xyz[j].cone, (k+2)*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
   xyz[j].cone[k] = atom2; xyz[j].cone[k+1] = 0;
   j=atom2-1; k=0; while(xyz[j].cone[k]!=0) k++;
   if((xyz[j].cone = (int *) realloc(xyz[j].cone, (k+2)*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
   xyz[j].cone[k] = atom1; xyz[j].cone[k+1] = 0;
   }
  }
 }
for(i=0;input[i]!=NULL;i++)
 {
 if((read_firststrcmp(input[i], "delbond"))!=NULL)
  {
  if((c = sscanf(input[i], "%*s%d%d", &atom1, &atom2))!=2) {puts("\nProblems when reading delbond data...\n"); exit(1);}
  if(atom1<1 || atom1>natoms) {printf("\nAtom number %d out of range for delbond...\n\n", atom1); exit(1);}
  if(atom2<1 || atom2>natoms) {printf("\nAtom number %d out of range for delbond...\n\n", atom2); exit(1);}
  j=atom1-1; k=0; while(xyz[j].cone[k]!=atom2 && xyz[j].cone[k]!=0) k++; while(xyz[j].cone[k]!=0) {xyz[j].cone[k] = xyz[j].cone[k+1]; k++;}
  j=atom2-1; k=0; while(xyz[j].cone[k]!=atom1 && xyz[j].cone[k]!=0) k++; while(xyz[j].cone[k]!=0) {xyz[j].cone[k] = xyz[j].cone[k+1]; k++;}
  }
 }
for(i=0;input[i]!=NULL;i++)
 {
 if((read_firststrcmp(input[i], "delallbonds"))!=NULL)
  {
  if((c = sscanf(input[i], "%*s%d", &atom1))!=1) {puts("\nProblems when reading delallbonds data...\n"); exit(1);}
  if(atom1<1 || atom1>natoms) {printf("\nAtom number %d out of range for delallbonds...\n\n", atom1); exit(1);}
  xyz[atom1-1].cone[0] = 0;
  for(j=0;j<natoms;j++) {k=0; while(xyz[j].cone[k]!=atom1 && xyz[j].cone[k]!=0) k++; while(xyz[j].cone[k]!=0) {xyz[j].cone[k] = xyz[j].cone[k+1]; k++;}}
  }
 }
}

/*********************************************************************************************************************************************************************************/
// functions for gauges
// reading gauges and making their settings (dashed is default style)
void atoms_make_gauge_settings(GAUGE *gauge, int ngauges, ATOMS *xyz, int natoms, int nmolecules, double *bondradius, char **input)
{
int i, j, k, c;

for(i=0;i<ngauges;i++)
 {
 k=0; for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "addgauge"))!=NULL) {if(k==i) break; else k++;}
 if((c = sscanf(input[j], "%*s%d%d", &gauge[i].a1, &gauge[i].a2))!=2) {puts("\nProblems when reading addgauge data...\n"); exit(1);}
 if(gauge[i].a1==gauge[i].a2) {puts("\nAtom numbers for making meaningful gauge must differ...\n"); exit(1);}
 if(gauge[i].a1<1 || gauge[i].a1>natoms) {printf("\nAtom number %d out of range for addgauge...\n\n", gauge[i].a1); exit(1);}
 if(gauge[i].a2<1 || gauge[i].a2>natoms) {printf("\nAtom number %d out of range for addgauge...\n\n", gauge[i].a2); exit(1);}
 if((c = sscanf(input[j], "%*s%*d%*d%lf", &gauge[i].scale))!=1) gauge[i].scale = 0.4;
 if((c = sscanf(input[j], "%*s%*d%*d%*f%*s%lf%lf%lf", &gauge[i].rgb[0], &gauge[i].rgb[1], &gauge[i].rgb[2]))!=3) for(k=0;k<3;k++) gauge[i].rgb[k] = 0.0;
 gauge[i].style = 1;
 if((strstr(input[j], "dotted"))!=NULL) gauge[i].style = 2;
 if((strstr(input[j], "dotdashed"))!=NULL) gauge[i].style = 3;
 if((strstr(input[j], "solid"))!=NULL) gauge[i].style = 4;
 if(xyz[gauge[i].a1-1].molid==xyz[gauge[i].a2-1].molid) gauge[i].rad = gauge[i].scale*bondradius[xyz[gauge[i].a1-1].ordid];
  else gauge[i].rad = gauge[i].scale*bondradius[nmolecules];
 }
}

