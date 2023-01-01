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
#include "charvib.h"
#include "charread.h"
#include "charcheck.h"
#include "charwrite.h"

/*********************************************************************************************************************************************************************************/
// functions for vibrations
// processing options from a line
void vib_process_options(int natoms, VIBRATION *vib, int num, char *line, int global)
{
int i, j, k, c, tookmore, a1, a2, ats[natoms+1];
char *string, *token, help[50];
double dval;

const char keywords[] = " vib vibrationscale arradscale conradscale atomiccolors vibrationcolor includeatoms excludeatoms ";

string = line;
token = read_next_token(string);
string = NULL;
while(token!=NULL)
 {
 strcpy(help, " "); strcat(help, token); strcat(help, " ");
 if((strstr(keywords, help))==NULL) {printf("\nUnknown keyword '%s' in vibration section...\n\n", token); exit(1);}
 tookmore = 0;
 if((strcmp(token, "vib"))==0 && global!=1)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  if((c = sscanf(token, "%d", &vib[num].number))!=1) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  }
 if((strcmp(token, "vibrationscale"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  if((c = sscanf(token, "%lf", &dval))!=1) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  if(global==1) for(i=0;i<num;i++) vib[i].scale *= dval;
   else vib[num].scale *= dval;
  }
 if((strcmp(token, "arradscale"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  if((c = sscanf(token, "%lf", &dval))!=1) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  if(global==1) for(i=0;i<num;i++) vib[i].arradscale = dval;
   else vib[num].arradscale = dval;
  }
 if((strcmp(token, "conradscale"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  if((c = sscanf(token, "%lf", &dval))!=1) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  if(global==1) for(i=0;i<num;i++) vib[i].conradscale = dval;
   else vib[num].conradscale = dval;
  }
 if((strcmp(token, "atomiccolors"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  if((strcmp(token, "on"))!=0 && (strcmp(token, "off"))!=0) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
  if((strcmp(token, "on"))==0)
   {
   if(global==1) for(i=0;i<num;i++) vib[i].atomiccolor = 1;
    else vib[num].atomiccolor = 1;
   }
  }
 if((strcmp(token, "vibrationcolor"))==0)
  {
  for(j=0;j<3;j++)
   {
   token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
   if((c = sscanf(token, "%lf", &dval))!=1) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
   if(global==1) for(i=0;i<num;i++) {vib[i].rgb[j] = dval; vib[i].atomiccolor = 0;}
    else {vib[num].rgb[j] = dval; vib[num].atomiccolor = 0;}
   }
  }
 if(token!=NULL && (strcmp(token, "includeatoms"))==0)
  {
  tookmore = 1;
  ats[0] = 0;
  while((token = read_next_token(string))!=NULL)
   {
   if((c = sscanf(token, "%d-%d", &a1, &a2))==2)
    {
    if(a1<1 || a1>natoms) {printf("\nAtom number %d out of range for vibration...\n\n", a1); exit(1);}
    if(a2<1 || a2>natoms) {printf("\nAtom number %d out of range for vibration...\n\n", a2); exit(1);}
    if(a1>a2) {i = a1; a1 = a2; a2 = i;}
    i = 0; while(ats[i]!=0) i++;
    while(a1<=a2) {ats[i] = a1; ats[i+1] = 0; a1++; i++;}
    }
   if(c==1)
    {
    if(a1<1 || a1>natoms) {printf("\nAtom number %d out of range for vibration...\n\n", a1); exit(1);}
    i = 0; while(ats[i]!=0) i++;
    ats[i] = a1; ats[i+1] = 0;
    }
   if(c==0) break;
   }
  if(global==1)
   {
   for(i=0;i<num;i++)
    {
    for(j=0;ats[j]!=0;j++)
     {
     k = 0; while(vib[i].atoms[k]!=0 && vib[i].atoms[k]!=ats[j]) k++;
     if(vib[i].atoms[k]==0) {vib[i].atoms[k] = ats[j]; vib[i].atoms[k+1] = 0;}
     }
    }
   }
   else
    {
    for(j=0;ats[j]!=0;j++)
     {
     k = 0; while(vib[num].atoms[k]!=0 && vib[num].atoms[k]!=ats[j]) k++;
     if(vib[num].atoms[k]==0) {vib[num].atoms[k] = ats[j]; vib[num].atoms[k+1] = 0;}
     }
    }
  }
 if(token!=NULL && (strcmp(token, "excludeatoms"))==0)
  {
  tookmore = 1;
  ats[0] = 0;
  while((token = read_next_token(string))!=NULL)
   {
   if((c = sscanf(token, "%d-%d", &a1, &a2))==2)
    {
    if(a1<1 || a1>natoms) {printf("\nAtom number %d out of range for vibration...\n\n", a1); exit(1);}
    if(a2<1 || a2>natoms) {printf("\nAtom number %d out of range for vibration...\n\n", a2); exit(1);}
    if(a1>a2) {i = a1; a1 = a2; a2 = i;}
    i = 0; while(ats[i]!=0) i++;
    while(a1<=a2) {ats[i] = a1; ats[i+1] = 0; a1++; i++;}
    }
   if(c==1)
    {
    if(a1<1 || a1>natoms) {printf("\nAtom number %d out of range for vibration...\n\n", a1); exit(1);}
    i = 0; while(ats[i]!=0) i++;
    ats[i] = a1; ats[i+1] = 0;
    }
   if(c==0) break;
   }
  if(global==1)
   {
   for(i=0;i<num;i++)
    {
    if(vib[i].atoms[0]==0) {for(j=0;j<natoms;j++) vib[i].atoms[j] = j+1; vib[i].atoms[j] = 0;}
    for(j=0;ats[j]!=0;j++)
     {
     k = 0; while(vib[i].atoms[k]!=0 && vib[i].atoms[k]!=ats[j]) k++;
     if(vib[i].atoms[k]==ats[j]) {while(vib[i].atoms[k]!=0) {vib[i].atoms[k] = vib[i].atoms[k+1]; k++;}}
     }
    }
   }
   else
    {
    if(vib[num].atoms[0]==0) {for(j=0;j<natoms;j++) vib[num].atoms[j] = j+1; vib[num].atoms[j] = 0;}
    for(j=0;ats[j]!=0;j++)
     {
     k = 0; while(vib[num].atoms[k]!=0 && vib[num].atoms[k]!=ats[j]) k++;
     if(vib[num].atoms[k]==ats[j]) {while(vib[num].atoms[k]!=0) {vib[num].atoms[k] = vib[num].atoms[k+1]; k++;}}
     }
    }
  }
 if(token!=NULL && tookmore!=1) token = read_next_token(string);
 }
}

// making arrows for a vibration
int vib_make_arrows(VIBRATION actvib, TYCKA *vibtyc, KONUS *vibkon, ATOMS *xyz, int natoms, double *bondradius)
{
int i, j, enlarge;
double dist;

enlarge = 0;

for(i=0;i<natoms;i++)
 {
 vibtyc[i].rad = 0.8*bondradius[vibtyc[i].ordid]*actvib.arradscale;
 vibkon[i].rad = 2.0*vibtyc[i].rad*actvib.conradscale;
 for(j=0;j<3;j++) vibkon[i].end[j] *= actvib.scale;
 dist = 0.0; for(j=0;j<3;j++) dist += vibkon[i].end[j]*vibkon[i].end[j]; dist = sqrt(dist);
 for(j=0;j<3;j++) vibtyc[i].end[j] = vibkon[i].end[j]*(1.0 - vibkon[i].rad/dist);
 for(j=0;j<3;j++) {vibtyc[i].end[j] += vibtyc[i].start[j]; vibkon[i].end[j] += vibtyc[i].start[j];}
 if(actvib.atomiccolor==1) {for(j=0;j<3;j++) vibtyc[i].rgb[j] = xyz[i].rgb[j];}
  else {for(j=0;j<3;j++) vibtyc[i].rgb[j] = actvib.rgb[j];}
 for(j=0;actvib.atoms[j]!=0;j++) if(actvib.atoms[j] == i+1 && dist > 0.0 && dist < vibkon[i].rad) enlarge = 1;
 }

return enlarge;
}

/*********************************************************************************************************************************************************************************/
// main function for vibrations
// making vibrations
void vib_make_vibrations(ATOMS *xyz, int natoms, int nmolecules, int centerofmass, double *centerofmasscoor, double convertxyzval,
                         int hideatoms, int atomcolor, int bondcolor, double *atomradiusscale, double *lightatomradiusscale,
                         double *bondradius, double **defatomcolor, double **defbondcolor, TYCKA *tyc, KONUS *kon, int narrows,
                         GAUGE *gauge, int ngauges, int camspec, int useangle, double *camera, int scenerot, double scenerotangle,
                         int camrot, double *camerarot, int camrq, double *camerarq, int camsky, double *camerasky, int camzoom,
                         double camerazoom, int handedness, int finish, double *finset, char **input, char **extinput,
                         int extfileformat, char *outfilename, char *suffix, int outformat, char *command, char *comend)
{
int i, j, k, l, nvibs, enlarge;
char actfilename[500], ordstr[20], actcommand[550];
FILE *outfile;
VIBRATION *vib;
TYCKA *vibtyc;
KONUS *vibkon;

// initialization
vib = NULL; vibtyc = NULL; vibkon = NULL;

// determining number of vibrations and allocating memory
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "vibration"))!=NULL) break; i++;
for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endvibration"))!=NULL) break;
nvibs = 0;
for(k=i;k<j;k++) if((strstr(input[k], "vib "))!=NULL) nvibs++;
if(nvibs==0) {puts("\nNo vibration specification found in 'vibration' section...\n"); exit(1);}
if((vib = (VIBRATION *) malloc(nvibs*sizeof(VIBRATION)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
for(k=0;k<nvibs;k++) if((vib[k].atoms = (int *) malloc((natoms+1)*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((vibtyc = (TYCKA *) malloc(natoms*sizeof(TYCKA)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((vibkon = (KONUS *) malloc(natoms*sizeof(KONUS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}

/*********************************************************************************************************************************************************************************/
// VARIOUS SETTINGS IMPORTANT TO USERS

// default global settings for vibrations
for(k=0;k<nvibs;k++)
 {
 vib[k].scale = 10.0;
 vib[k].atomiccolor = 0;
 for(l=0;l<3;l++) vib[k].rgb[l] = 0.0;
 vib[k].arradscale = 1.0;
 vib[k].conradscale = 1.0;
 vib[k].atoms[0] = 0;
 }

// modification of global settings for vibrations specified by user
for(k=i;k<j;k++)
 {
 if((strstr(input[k], "vib "))==NULL) vib_process_options(natoms, vib, nvibs, input[k], 1);
 }

// modification of settings for particular vibrations specified by user
// order number of vibration is filled during this
l = 0;
for(k=i;k<j;k++)
 {
 if((strstr(input[k], "vib "))!=NULL)
  {
  vib_process_options(natoms, vib, l, input[k], 0);
  l++;
  }
 }

// finishing the settings
for(i=0;i<nvibs;i++)
 {
 if(vib[i].atoms[0]==0)
  {
  for(j=0;j<natoms;j++) vib[i].atoms[j] = j+1;
  vib[i].atoms[j] = 0;
  }
 }
/*********************************************************************************************************************************************************************************/

// reading vibrations, making the arrows, writing output files and possibly running povray in a loop
for(i=0;i<nvibs;i++)
 {
 // reading vibration and making arrows
 read_vibration_molden(vibtyc, vibkon, xyz, natoms, centerofmass, centerofmasscoor, vib[i].number, extinput);
 enlarge = vib_make_arrows(vib[i], vibtyc, vibkon, xyz, natoms, bondradius);
 // checking vibration arrows
 for(j=0;vib[i].atoms[j]!=0;j++)
  {
  check_single_arrow(vibtyc[vib[i].atoms[j]-1], vibkon[vib[i].atoms[j]-1], vib[i].atoms[j], vib[i].number);
  }
 if(enlarge == 1) printf("Enlarge scale of vibration %d or exclude atoms with small movements in order to display all arrows correctly...\n", vib[i].number);
 // writing output
 strcpy(actfilename, outfilename);
 strcat(actfilename, ".vib.");
 sprintf(ordstr, "%d", vib[i].number);
 strcat(actfilename, ordstr);
 strcat(actfilename, suffix);
 if((outfile = fopen(actfilename, "w"))==NULL) {printf("\nOpening file %s failed!\n\n", actfilename); exit(1);}
 // common part of output
 outfile = write_basic_settings(xyz, natoms, tyc, kon, narrows, camspec, useangle, camera, convertxyzval, scenerot,
                                scenerotangle, camrot, camerarot, camrq, camerarq, camsky, camerasky, camzoom, camerazoom,
                                handedness, finish, finset, outfile, outformat);
 outfile = write_xyz(xyz, natoms, nmolecules, hideatoms, atomcolor, bondcolor, finish, atomradiusscale,
                     lightatomradiusscale, bondradius, defatomcolor, defbondcolor, outfile, outformat);
 outfile = write_gauges(xyz, natoms, gauge, ngauges, hideatoms, finish, atomradiusscale,
                        lightatomradiusscale, bondradius, outfile, outformat);
 outfile = write_arrows(tyc, kon, narrows, finish, outfile, outformat);
 if(narrows!=0) puts("Picture of vibration contains also arrows specified in charmol input...");
 // writing vibration
 for(j=0;vib[i].atoms[j]!=0;j++)
  {
  outfile = write_single_arrow(vibtyc[vib[i].atoms[j]-1], vibkon[vib[i].atoms[j]-1], finish, outfile, outformat);
  }
 if((fclose(outfile))==EOF) {printf("\nClosing file %s failed!\n\n", actfilename); exit(1);}
 // running povray
 if(outformat==1 && (read_firststrcmp(command, "norendering"))==NULL)
  {
  strcpy(actcommand, command);
  strcat(actcommand, actfilename);
  strcat(actcommand, comend);
  write_run_system_command(actcommand);
  }
 }

// deallocate memory
for(i=0;i<nvibs;i++) free(vib[i].atoms); free(vib); vib = NULL;
free(vibtyc); vibtyc = NULL;
free(vibkon); vibkon = NULL;
}

