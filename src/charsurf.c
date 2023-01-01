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
#include "charsurf.h"
#include "charread.h"
#include "charcheck.h"
#include "charwrite.h"
#include "charmesh.h"

/*********************************************************************************************************************************************************************************/
// functions for surfaces
// processing global options from a line
void surf_process_options(SURFACE *surf, int num, char *line, int global)
{
int i, j, c, tookmore, s1, s2, *surfs;
char *string, *token, keyhelp[50];
double dval;

const char keywords[] = " cubefile longfilenames isovalue color positivecolor negativecolor transparency colormapped style orbitals ";

surfs = NULL;
if((surfs = (int *) malloc(1*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
surfs[0] = 0;

string = line;
token = read_next_token(string);
string = NULL;
while(token!=NULL)
 {
 strcpy(keyhelp, " "); strcat(keyhelp, token); strcat(keyhelp, " ");
 if((strstr(keywords, keyhelp))==NULL) {printf("\nUnknown keyword '%s' in surface section...\n\n", token); exit(1);}
 tookmore = 0;
 if((strcmp(token, "cubefile"))==0 && global!=1)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  if((c = sscanf(token, "%s", surf[num].filename))!=1) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  }
 if((strcmp(token, "longfilenames"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  if((strcmp(token, "on"))!=0 && (strcmp(token, "off"))!=0) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  if((strcmp(token, "on"))==0)
   {
   if(global==1) for(i=0;i<num;i++) surf[i].longname = 1;
    else surf[num].longname = 1;
   }
   else
    {
    if(global==1) for(i=0;i<num;i++) surf[i].longname = 0;
     else surf[num].longname = 0;
    }
  }
 if((strcmp(token, "isovalue"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  if((c = sscanf(token, "%lf", &dval))!=1) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  if(global==1) for(i=0;i<num;i++) surf[i].isoval = dval;
   else surf[num].isoval = dval;
  }
 if((strcmp(token, "color"))==0)
  {
  for(j=0;j<3;j++)
   {
   token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
   if((c = sscanf(token, "%lf", &dval))!=1) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
   if(global==1) for(i=0;i<num;i++) surf[i].rgb[j] = dval;
    else surf[num].rgb[j] = dval;
   }
  }
 if((strcmp(token, "positivecolor"))==0)
  {
  for(j=0;j<3;j++)
   {
   token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
   if((c = sscanf(token, "%lf", &dval))!=1) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
   if(global==1) for(i=0;i<num;i++) surf[i].rgbpos[j] = dval;
    else surf[num].rgbpos[j] = dval;
   }
  }
 if((strcmp(token, "negativecolor"))==0)
  {
  for(j=0;j<3;j++)
   {
   token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
   if((c = sscanf(token, "%lf", &dval))!=1) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
   if(global==1) for(i=0;i<num;i++) surf[i].rgbneg[j] = dval;
    else surf[num].rgbneg[j] = dval;
   }
  }
 if((strcmp(token, "transparency"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  if((c = sscanf(token, "%lf", &dval))!=1) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  if(global==1) for(i=0;i<num;i++) surf[i].transparency = dval;
   else surf[num].transparency = dval;
  }
 if((strcmp(token, "colormapped"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  if(global==1)
   {
   for(i=0;i<num;i++)
    {
    surf[i].colmap = 1;
    if((c = sscanf(token, "%s", surf[i].potfilename))!=1) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
    }
   }
   else
    {
    surf[num].colmap = 1;
    if((c = sscanf(token, "%s", surf[num].potfilename))!=1) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
    }
  }
 if((strcmp(token, "style"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  c = 0;
  if((strcmp(token, "solid"))==0)      c = 1;
  if((strcmp(token, "wireframe"))==0)  c = 2;
  if(c==0) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  if(global==1) for(i=0;i<num;i++) surf[i].style = c;
   else surf[num].style = c;
  }
 if((strcmp(token, "orbitals"))==0)
  {
  tookmore = 1;
  while((token = read_next_token(string))!=NULL)
   {
   if((c = sscanf(token, "%d-%d", &s1, &s2))==2)
    {
    if(s1>s2) {i = s1; s1 = s2; s2 = i;}
    j = 0; while(surfs[j]!=0) j++;
    if((surfs = (int *) realloc(surfs, ((s2-s1+1)+j+1)*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
    for(i=s1;i<=s2;i++,j++) surfs[j] = i; surfs[j] = 0;
    }
   if(c==1)
    {
    j = 0; while(surfs[j]!=0) j++;
    if((surfs = (int *) realloc(surfs, (1+j+1)*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
    surfs[j] = s1; surfs[j+1] = 0;
    }
   if(c==0) break;
   }
  c = 0; while(surfs[c]!=0) c++;
  if(c==0) {puts("\nNo orbitals specified in surface section...\n"); exit(1);}
  if(global==1)
   {
   for(i=0;i<num;i++)
    {
    if((surf[i].surfs = (int *) realloc(surf[i].surfs, c*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
    for(j=0;j<c;j++) surf[i].surfs[j] = surfs[j];
    surf[i].nsurf = c;
    }
   }
   else
    {
    if((surf[num].surfs = (int *) realloc(surf[num].surfs, c*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
    for(j=0;j<c;j++) surf[num].surfs[j] = surfs[j];
    surf[num].nsurf = c;
    }
  }
 if(token!=NULL && tookmore!=1) token = read_next_token(string);
 }

free(surfs); surfs = NULL;
}

/*********************************************************************************************************************************************************************************/
// main function for surfaces
// making surfaces
void surf_make_surfaces(ATOMS *xyz, int natoms, int nmolecules, int centerofmass, double *centerofmasscoor, double convertxyzval,
                        int hideatoms, int atomcolor, int bondcolor, double *atomradiusscale, double *lightatomradiusscale,
                        double *bondradius, double **defatomcolor, double **defbondcolor, TYCKA *tyc, KONUS *kon, int narrows,
                        GAUGE *gauge, int ngauges, int camspec, int useangle, double *camera, int scenerot, double scenerotangle,
                        int camrot, double *camerarot, int camrq, double *camerarq, int camsky, double *camerasky, int camzoom,
                        double camerazoom, int handedness, int finish, double *finset, char **input, char **extinput,
                        int extfileformat, char *outfilename, char *suffix, int outformat, char *command, char *comend)
{
int i, j, k, l, s, nsurfs, *gpoints, totpoints, totsurfs;
char **surfinput, **potinput, actfilename[500], actcommand[550], namehelp[300];
double *gmin, *griddist, *gridx, *gridy, *gridz, **gridval, *gridpot;
FILE *outfile;
SURFACE *surf;
MESH mesh;

// initialization
surf = NULL; surfinput = NULL; potinput = NULL; gpoints = NULL; gmin = NULL; griddist = NULL;
gridx = NULL; gridy = NULL; gridz = NULL; gridval = NULL; gridpot = NULL;

// determining number of surfaces
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "surface"))!=NULL) break; i++;
for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endsurface"))!=NULL) break;
nsurfs = 0;
for(k=i;k<j;k++) if((strstr(input[k], "cubefile"))!=NULL) nsurfs++;
if(nsurfs==0) {puts("\nNo cubefile specification found in 'surface' section...\n"); exit(1);}
if((surf = (SURFACE *) malloc(nsurfs*sizeof(SURFACE)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((gpoints = (int *) malloc(3*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((gmin = (double *) malloc(3*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((griddist = (double *) malloc(3*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}

/*********************************************************************************************************************************************************************************/
// VARIOUS SETTINGS IMPORTANT TO USERS

// default settings
for(k=0;k<nsurfs;k++)
 {
 surf[k].isoval = 0.0;
 surf[k].sign = 0;
 surf[k].transparency = 0.0;
 surf[k].rgb[0] = 1.0; surf[k].rgb[1] = 1.0; surf[k].rgb[2] = 0.0;
 surf[k].rgbpos[0] = 0.0; surf[k].rgbpos[1] = 0.0; surf[k].rgbpos[2] = 1.0;
 surf[k].rgbneg[0] = 1.0; surf[k].rgbneg[1] = 0.0; surf[k].rgbneg[2] = 0.0;
 surf[k].nsurf = 0;
 surf[k].surfs = NULL;
 surf[k].norb = 0;
 surf[k].orbs = NULL;
 surf[k].longname = 0;
 surf[k].colmap = 0;
 surf[k].style = 1;
 }

// modification of global settings for surfaces specified by user
for(k=i;k<j;k++)
 {
 if((strstr(input[k], "cubefile"))==NULL) surf_process_options(surf, nsurfs, input[k], 1);
 }

// modification of settings for particular surfaces specified by user
// cube filenames are filled during this
l = 0;
for(k=i;k<j;k++)
 {
 if((strstr(input[k], "cubefile"))!=NULL)
  {
  surf_process_options(surf, l, input[k], 0);
  l++;
  }
 }
/*********************************************************************************************************************************************************************************/

// checking settings of surfaces for hard problems (checking all because sign is not known yet)
for(s=0;s<nsurfs;s++)
 {
 check_mesh(surf[s].rgb, surf[s].transparency, surf[s].isoval, 1);
 check_mesh(surf[s].rgbpos, surf[s].transparency, surf[s].isoval, 1);
 check_mesh(surf[s].rgbneg, surf[s].transparency, -surf[s].isoval, 1);
 }

// reading and depiction of all required surfaces in a loop
for(s=0;s<nsurfs;s++)
 {
 // reading cube file
 surfinput = read_file_to_char(surf[s].filename);
 // checking geometry, required surfaces and finishing their settings
 // if user did not specify orbitals (surfaces) for particular cube file, all of them will be done
 read_check_geom_surf_cube(surf, s, xyz, natoms, centerofmass, centerofmasscoor, surfinput);
 if(surf[s].sign==1 && fabs(surf[s].isoval)<1E-10) surf[s].isoval = 0.1;
 if(surf[s].sign==2 && fabs(surf[s].isoval)<1E-10) surf[s].isoval = 0.05;
 if(surf[s].norb==0) totsurfs = surf[s].nsurf;
  else totsurfs = surf[s].norb;
 // reading grid information, preparing grid and reading volumetric data
 read_grid_info_cube(gpoints, gmin, griddist, centerofmass, centerofmasscoor, surfinput);
 totpoints = gpoints[0]*gpoints[1]*gpoints[2];
 if((gridx = (double *) realloc(gridx, totpoints*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((gridy = (double *) realloc(gridy, totpoints*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((gridz = (double *) realloc(gridz, totpoints*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((gridval = (double **) malloc(totsurfs*sizeof(double *)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 for(i=0;i<totsurfs;i++) if((gridval[i] = (double *) malloc(totpoints*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
#ifdef _OPENMP
 #pragma omp parallel for default(shared) private(j, k, l)
#endif
 for(i=0;i<gpoints[0];i++)
  {
  for(j=0;j<gpoints[1];j++)
   {
   for(k=0;k<gpoints[2];k++)
    {
    l = i*gpoints[1]*gpoints[2] + j*gpoints[2] + k;
    gridx[l] = gmin[0] + i*griddist[0];
    gridy[l] = gmin[1] + j*griddist[1];
    gridz[l] = gmin[2] + k*griddist[2];
    }
   }
  }
 read_grid_values_cube(gridval, totpoints, totsurfs, surfinput);
 // checking signs of actual grid values (difference density - two signs although it is density, not orbital)
 // code could be possibly modified to use two signs for everything (orbitals and also all densities etc.)
 read_check_grid_value_signs_cube(surf, s, gridval, totpoints, totsurfs);
 // checking things for colormapped surface and reading potential
 if(surf[s].colmap==1)
  {
  potinput = read_file_to_char(surf[s].potfilename);
  read_check_geom_grid_cube(xyz, natoms, centerofmass, centerofmasscoor, gpoints, gmin, griddist, potinput);
  if((gridpot = (double *) realloc(gridpot, totpoints*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
  read_grid_values_oneset_cube(gridpot, totpoints, potinput);
  }
 // making required surfaces from one cube file
 for(i=0;i<surf[s].nsurf;i++)
  {
  // writing output
  strcpy(actfilename, surf[s].filename);
  j = 0; while(actfilename[j]!='\0') j++;
  if(actfilename[j-1]=='e' && actfilename[j-2]=='b' && actfilename[j-3]=='u' && actfilename[j-4]=='c' && actfilename[j-5]=='.') actfilename[j-5] = '\0';
  if(surf[s].norb!=0) {sprintf(namehelp, "%d", surf[s].orbs[surf[s].surfs[i]]); strcat(actfilename, ".orb."); strcat(actfilename, namehelp);}
  if(surf[s].longname==1)
   {
   if(surf[s].sign==1) sprintf(namehelp, ".iso-%f.rgb-%.3f-%.3f-%.3f.transp-%.3f",
                               surf[s].isoval, surf[s].rgb[0], surf[s].rgb[1], surf[s].rgb[2], surf[s].transparency);
    else sprintf(namehelp, ".iso-%f.posrgb-%.3f-%.3f-%.3f.negrgb-%.3f-%.3f-%.3f.transp-%.3f",
                 surf[s].isoval, surf[s].rgbpos[0], surf[s].rgbpos[1], surf[s].rgbpos[2],
                 surf[s].rgbneg[0], surf[s].rgbneg[1], surf[s].rgbneg[2], surf[s].transparency);
   strcat(actfilename, namehelp);
   }
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
  // generating and writing mesh
  switch(surf[s].sign)
   {
   case 1:
    if(surf[s].colmap==1) mesh = mesh_generate_mesh(gpoints, griddist, gridx, gridy, gridz, gridval[surf[s].surfs[i]], gridpot, surf[s].isoval);
     else mesh = mesh_generate_mesh(gpoints, griddist, gridx, gridy, gridz, gridval[surf[s].surfs[i]], NULL, surf[s].isoval);
    switch(surf[s].style)
     {
     case 1: outfile = write_mesh_solid(mesh, surf[s].rgb, surf[s].transparency, surf[s].isoval, finish, outfile, outformat); break;
     case 2: outfile = write_mesh_wireframe(mesh, surf[s].rgb, surf[s].transparency, surf[s].isoval, finish, outfile, outformat); break;
     }
    mesh = mesh_clean_mesh(mesh);
    break;
   case 2:
    if(surf[s].colmap==1) mesh = mesh_generate_mesh(gpoints, griddist, gridx, gridy, gridz, gridval[surf[s].surfs[i]], gridpot, surf[s].isoval);
     else mesh = mesh_generate_mesh(gpoints, griddist, gridx, gridy, gridz, gridval[surf[s].surfs[i]], NULL, surf[s].isoval);
    switch(surf[s].style)
     {
     case 1: outfile = write_mesh_solid(mesh, surf[s].rgbpos, surf[s].transparency, surf[s].isoval, finish, outfile, outformat); break;
     case 2: outfile = write_mesh_wireframe(mesh, surf[s].rgbpos, surf[s].transparency, surf[s].isoval, finish, outfile, outformat); break;
     }
    mesh = mesh_clean_mesh(mesh);
    if(surf[s].colmap==1) mesh = mesh_generate_mesh(gpoints, griddist, gridx, gridy, gridz, gridval[surf[s].surfs[i]], gridpot, -surf[s].isoval);
     else mesh = mesh_generate_mesh(gpoints, griddist, gridx, gridy, gridz, gridval[surf[s].surfs[i]], NULL, -surf[s].isoval);
    switch(surf[s].style)
     {
     case 1: outfile = write_mesh_solid(mesh, surf[s].rgbneg, surf[s].transparency, -surf[s].isoval, finish, outfile, outformat); break;
     case 2: outfile = write_mesh_wireframe(mesh, surf[s].rgbneg, surf[s].transparency, -surf[s].isoval, finish, outfile, outformat); break;
     }
    mesh = mesh_clean_mesh(mesh);
    break;
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
 for(i=0;i<totsurfs;i++) free(gridval[i]); free(gridval); gridval = NULL;
 for(i=0;surfinput[i]!=NULL;i++) free(surfinput[i]); free(surfinput); surfinput = NULL;
 if(potinput!=NULL) {for(i=0;potinput[i]!=NULL;i++) free(potinput[i]); free(potinput); potinput = NULL;}
 }

// deallocate memory
free(surf); surf = NULL;
free(gpoints); gpoints = NULL;
free(gmin); gmin = NULL;
free(griddist); griddist = NULL;
free(gridx); gridx = NULL;
free(gridy); gridy = NULL;
free(gridz); gridz = NULL;
if(gridpot!=NULL) {free(gridpot); gridpot = NULL;}
}

