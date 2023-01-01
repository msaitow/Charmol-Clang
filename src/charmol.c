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
#include "charhelp.h"
#include "charread.h"
#include "charcheck.h"
#include "charwrite.h"
#include "charatoms.h"
#include "charatomsinfo.h"
#include "chararrow.h"
#include "charvib.h"
#include "charorb.h"
#include "charsurf.h"

/*********************************************************************************************************************************************************************************/
// main program of CHARMOL
int main(int argc, char **argv)
{

if(argc<2) {printf("Syntax: charmol inputfile [povrayverbose] [hiderendering/norendering]\nTry '-h' or '--help' option for help.\n"); exit(1);}
if((strcmp(argv[1], "--help"))==0 || (strcmp(argv[1], "-h"))==0) {help_print(); exit(0);}

int i, j, k, l, c, natoms, finish, narrows, atomcolor, bondcolor, hideatoms, camspec, camrot, camrq, camsky, camzoom, nmolecules, lastmolid, extfileformat, useangle;
int atomsinfo, centerofmass, handedness, moldenview, scenerot, help, ngauges, outformat, render, hide, verbose, *ids, vibrations, orbitals, surfaces;
char **input, **extinput, **atoms, atomsinfofilename[500], extfilename[500], handspec[10], outfilename[500], suffix[10], command[550], comend[20], *string, *token;
double camera[3], camerarot[3], camerarq[4], camerasky[3], camerazoom, scenerotangle, finset[5], val, vals[3], convertxyzval, centerofmasscoor[3], sumofmass;
double *bondradius, *arrowradius, *atomradiusscale, *lightatomradiusscale, *bondradiusscale, *arrowradiusscale, *coneradiusscale, **defatomcolor, **defbondcolor;
FILE *outfile;
ATOMS *xyz;
TYCKA *tyc;
KONUS *kon;
GAUGE *gauge;

// initialization
input = NULL; extinput = NULL; atoms = NULL; xyz = NULL; tyc = NULL; kon = NULL; gauge = NULL;
bondradius = NULL; arrowradius = NULL; atomradiusscale = NULL; lightatomradiusscale = NULL; bondradiusscale = NULL;
arrowradiusscale = NULL; coneradiusscale = NULL; defatomcolor = NULL; defbondcolor = NULL; ids = NULL;

// reading input file to input char matrix without comment lines (beginning with #)
input = read_file_to_char_withoutcomments(argv[1]);

// preparing basename for outputs
strcpy(outfilename, argv[1]);
i = 0; while(outfilename[i]!='\0') i++;
if(outfilename[i-1]=='p')
 {i -= 2; if(outfilename[i-2]=='.' && outfilename[i-1]=='i' && outfilename[i]=='n') outfilename[i-2] = '\0';}
if(outfilename[i-1]=='z')
 {i -= 2; if(outfilename[i-2]=='.' && outfilename[i-1]=='x' && outfilename[i]=='y') outfilename[i-2] = '\0';}

// reading external file, if specified
// external file format: 0 - no external file, 1 - molden, 2 - fchk, 3 - cube
extfileformat = 0;
for(i=0;input[i]!=NULL;i++)
 {
 if((read_firststrcmp(input[i], "moldenfile"))!=NULL)
  {
  extfileformat = 1;
  if((c = sscanf(input[i], "%*s%s", extfilename))!=1) {puts("\nProblems when reading molden filename...\n"); exit(1);}
  break;
  }
 if((read_firststrcmp(input[i], "fchkfile"))!=NULL)
  {
  extfileformat = 2;
  if((c = sscanf(input[i], "%*s%s", extfilename))!=1) {puts("\nProblems when reading fchk filename...\n"); exit(1);}
  break;
  }
 if((strstr(input[i], "cubefile"))!=NULL)
  {
  extfileformat = 3;
  string = &input[i][0]; token = read_next_token(string); string = NULL;
  while(token!=NULL && (strcmp(token, "cubefile"))!=0) token = read_next_token(string);
  if(token==NULL) {puts("\nProblems when reading cubefile specification...\n"); exit(1);}
  token = read_next_token(string);
  if(token==NULL || (c = sscanf(token, "%s", extfilename))!=1) {puts("\nProblems when reading cube filename...\n"); exit(1);}
  while(token!=NULL) token = read_next_token(string);
  break;
  }
 }
if(extfileformat!=0)
 {
 if(extfileformat==3) extinput = read_xyzpartofcubefile_to_char(extfilename);
  else extinput = read_file_to_char(extfilename);
 }

// reading xyz, setting number of atoms and molecules, rescaling xyz
switch(extfileformat)
 {
 case 0:
  xyz = read_xyz_input(input);
  break;
 case 1:
  xyz = read_xyz_molden(extinput);
  break;
 case 2:
  xyz = read_xyz_fchk(extinput);
  break;
 case 3:
  xyz = read_xyz_cube(extinput);
  break;
 }
natoms = 0; while(xyz[natoms].atomnum!=-1) natoms++;
nmolecules = 1; lastmolid = xyz[0].molid;
for(i=1;i<natoms;i++)
 {
 if(xyz[i].molid!=lastmolid)
  {
  help = 0;
  for(j=0;j<i;j++) if(xyz[j].molid==xyz[i].molid) help = 1;
  if(help==0) nmolecules++;
  lastmolid = xyz[i].molid;
  }
 }
convertxyzval = atoms_check_xyz_units(extfileformat, input, extinput);
atoms_rescale_xyz(xyz, natoms, convertxyzval);

// saving ids in user specified order and preparing data arrays for stuff, which can differ for different molecules
// on last position of arrays (for example bondradius[nmolecules]) the default will be stored for further use
if((ids = (int *) malloc(nmolecules*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((bondradius = (double *) malloc((nmolecules+1)*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((arrowradius = (double *) malloc((nmolecules+1)*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((atomradiusscale = (double *) malloc((nmolecules+1)*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((lightatomradiusscale = (double *) malloc((nmolecules+1)*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((bondradiusscale = (double *) malloc((nmolecules+1)*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((arrowradiusscale = (double *) malloc((nmolecules+1)*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((coneradiusscale = (double *) malloc((nmolecules+1)*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((defatomcolor = (double **) malloc((nmolecules+1)*sizeof(double *)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if((defbondcolor = (double **) malloc((nmolecules+1)*sizeof(double *)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
for(i=0;i<=nmolecules;i++)
 {
 if((defatomcolor[i] = (double *) malloc(3*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((defbondcolor[i] = (double *) malloc(3*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 }
i = 0; j = 0; ids[j] = xyz[i].molid; xyz[i].ordid = j; j++;
for(i=1;i<natoms;i++)
 {
 k=0; while(k<j && ids[k]!=xyz[i].molid) k++;
 if(k==j) {ids[j] = xyz[i].molid; xyz[i].ordid = j; j++;}
  else xyz[i].ordid = k;
 }

/*********************************************************************************************************************************************************************************/
// VARIOUS SETTINGS IMPORTANT TO USERS

// default settings
render = 1;
hide = 0;
verbose = 0;
outformat = 1;
atomcolor = 1;
bondcolor = 1;
hideatoms = 0;
atomsinfo = 0;
centerofmass = 1;
camspec = 0;
camrot = 0;
camrq = 0;
camsky = 0;
camzoom = 0;
useangle = 1;
handedness = -1;
scenerot = 0;
moldenview = 0;
finish = 1; finset[0] = 0.2; finset[1] = 0.8; finset[2] = 0.8; finset[3] = 0.2; finset[4] = 0.25;
for(i=0;i<=nmolecules;i++)
 {
 atomradiusscale[i] = 0.4;
 lightatomradiusscale[i] = 1.5;
 bondradiusscale[i] = 1.0;
 arrowradiusscale[i] = 1.0;
 coneradiusscale[i] = 1.0;
 defatomcolor[i][0] = 0.0; defatomcolor[i][1] = 0.0; defatomcolor[i][2] = 0.0;
 defbondcolor[i][0] = 0.7; defbondcolor[i][1] = 0.7; defbondcolor[i][2] = 0.7;
 }

// modification of settings specified by user
if(argc>2)
 {
 if((strcmp(argv[2], "norendering"))==0) render = 0;
 if((strcmp(argv[2], "hiderendering"))==0) hide = 1;
 if((strcmp(argv[2], "povrayverbose"))==0) verbose = 1;
 if(argc>3)
  {
  if((strcmp(argv[3], "norendering"))==0) render = 0;
  if((strcmp(argv[3], "hiderendering"))==0) hide = 1;
  if((strcmp(argv[3], "povrayverbose"))==0) verbose = 1;
  }
 }
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "atomradiusscale"))!=NULL)
 {if((c = sscanf(input[i], "%*s%lf%d", &val, &j))<1) {puts("\nProblems when reading atomradiusscale data...\n"); exit(1);}
  if(c==1) {for(k=0;k<=nmolecules;k++) atomradiusscale[k] = val; break;} else {k=0; while(ids[k]!=j && k<nmolecules) k++; if(k<nmolecules) atomradiusscale[k] = val;}}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "lightatomradiusscale"))!=NULL)
 {if((c = sscanf(input[i], "%*s%lf%d", &val, &j))<1) {puts("\nProblems when reading lightatomradiusscale data...\n"); exit(1);}
  if(c==1) {for(k=0;k<=nmolecules;k++) lightatomradiusscale[k] = val; break;} else {k=0; while(ids[k]!=j && k<nmolecules) k++; if(k<nmolecules) lightatomradiusscale[k] = val;}}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "bondradiusscale"))!=NULL)
 {if((c = sscanf(input[i], "%*s%lf%d", &val, &j))<1) {puts("\nProblems when reading bondradiusscale data...\n"); exit(1);}
  if(c==1) {for(k=0;k<=nmolecules;k++) bondradiusscale[k] = val; break;} else {k=0; while(ids[k]!=j && k<nmolecules) k++; if(k<nmolecules) bondradiusscale[k] = val;}}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "arrowradiusscale"))!=NULL)
 {if((c = sscanf(input[i], "%*s%lf%d", &val, &j))<1) {puts("\nProblems when reading arrowradiusscale data...\n"); exit(1);}
  if(c==1) {for(k=0;k<=nmolecules;k++) arrowradiusscale[k] = val; break;} else {k=0; while(ids[k]!=j && k<nmolecules) k++; if(k<nmolecules) arrowradiusscale[k] = val;}}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "coneradiusscale"))!=NULL)
 {if((c = sscanf(input[i], "%*s%lf%d", &val, &j))<1) {puts("\nProblems when reading coneradiusscale data...\n"); exit(1);}
  if(c==1) {for(k=0;k<=nmolecules;k++) coneradiusscale[k] = val; break;} else {k=0; while(ids[k]!=j && k<nmolecules) k++; if(k<nmolecules) coneradiusscale[k] = val;}}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "defaultatomcolor"))!=NULL)
 {atomcolor = 0; if((c = sscanf(input[i], "%*s%lf%lf%lf%d", &vals[0], &vals[1], &vals[2], &j))<3) {puts("\nProblems when reading defaultatomcolor data...\n"); exit(1);}
  if(c==3) {for(k=0;k<=nmolecules;k++) for(l=0;l<3;l++) defatomcolor[k][l] = vals[l]; break;}
   else {k=0; while(ids[k]!=j && k<nmolecules) k++; if(k<nmolecules) for(l=0;l<3;l++) defatomcolor[k][l] = vals[l];}}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "defaultbondcolor"))!=NULL)
 {bondcolor = 0; if((c = sscanf(input[i], "%*s%lf%lf%lf%d", &vals[0], &vals[1], &vals[2], &j))<3) {puts("\nProblems when reading defaultbondcolor data...\n"); exit(1);}
  if(c==3) {for(k=0;k<=nmolecules;k++) for(l=0;l<3;l++) defbondcolor[k][l] = vals[l]; break;}
   else {k=0; while(ids[k]!=j && k<nmolecules) k++; if(k<nmolecules) for(l=0;l<3;l++) defbondcolor[k][l] = vals[l];}}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "cameraposition"))!=NULL)
 {camspec = 1; if((c = sscanf(input[i], "%*s%lf%lf%lf", &camera[0], &camera[1], &camera[2]))!=3) {puts("\nProblems when reading cameraposition data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "camerarotation"))!=NULL)
 {camrot = 1; if((c = sscanf(input[i], "%*s%lf%lf%lf", &camerarot[0], &camerarot[1], &camerarot[2]))!=3) {puts("\nProblems when reading camerarotation data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "camerarotquat"))!=NULL)
 {camrq = 1; if((c = sscanf(input[i], "%*s%lf%lf%lf%lf", &camerarq[1], &camerarq[2], &camerarq[3], &camerarq[0]))!=4) {puts("\nProblems when reading camerarotquat data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "camerasky"))!=NULL)
 {camsky = 1; if((c = sscanf(input[i], "%*s%lf%lf%lf", &camerasky[0], &camerasky[1], &camerasky[2]))!=3) {puts("\nProblems when reading camerasky data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "camerazoom"))!=NULL)
 {camzoom = 1; if((c = sscanf(input[i], "%*s%lf", &camerazoom))!=1) {puts("\nProblems when reading camerazoom data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "finish"))!=NULL)
 {if((c = sscanf(input[i], "%*s%lf%lf%lf", &finset[0], &finset[1], &finset[2]))!=3) {puts("\nProblems when reading finish data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "ambient"))!=NULL)
 {if((c = sscanf(input[i], "%*s%lf", &finset[3]))!=1) {puts("\nProblems when reading ambient data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "shininess"))!=NULL)
 {if((c = sscanf(input[i], "%*s%lf", &finset[4]))!=1) {puts("\nProblems when reading shininess data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "atomsinfofile"))!=NULL)
 {atomsinfo = 1; if((c = sscanf(input[i], "%*s%s", atomsinfofilename))!=1) {puts("\nProblems when reading atomsinfofile filename...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "handedness"))!=NULL)
 {handedness = 0; if((c = sscanf(input[i], "%*s%s", handspec))!=1) {puts("\nProblems when reading handedness data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "scenerotation"))!=NULL)
 {scenerot = 1; if((c = sscanf(input[i], "%*s%lf", &scenerotangle))!=1) {puts("\nProblems when reading scenerotation data...\n"); exit(1);} break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "nofinish"))!=NULL) {finish = 0; break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "noatomcolor"))!=NULL) {atomcolor = 0; break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "nobondcolor"))!=NULL) {bondcolor = 0; break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "hideatoms"))!=NULL) {hideatoms = 1; break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "noangleuse"))!=NULL) {useangle = 0; break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "nocenterofmass"))!=NULL) {centerofmass = 0; break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "moldenview"))!=NULL) {moldenview = 1; break;}
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "output"))!=NULL) if((strstr(input[i], "vrml"))!=NULL) outformat = 2;

// finishing the settings
for(i=0;i<=nmolecules;i++)
 {
 bondradius[i] = 0.31*angstrom_to_bohr*atomradiusscale[i]*bondradiusscale[i];
 arrowradius[i] = 0.8*bondradius[i]*arrowradiusscale[i];
 }
if(atomsinfo==0) atoms = atomsinfo_fill_default_settings();
 else atoms = read_file_to_char_withoutcomments(atomsinfofilename);
if(handedness==0) {if((strstr(handspec, "left"))!=NULL) handedness = 1; else handedness = -1;}
if(moldenview==1) {handedness = -1; scenerot = 1; scenerotangle = -90.0; camspec = 0; camrot = 0;}
if(camrq==1) {scenerot = 0; camspec = 0; camrot = 0; camsky = 0;}
if(camsky==1) scenerot = 0;
if(camzoom==1) camerazoom = 1.0 - camerazoom/100.0;
if(outformat==1)
 {
 strcpy(suffix, ".pov");
 if(render==1) {if(hide==0) strcpy(command, "povray "); else strcpy(command, "povray -D ");}
  else strcpy(command, "norendering ");
 if(verbose==0) strcpy(comend, " >& /dev/null"); else strcpy(comend, "");
 }
 else strcpy(suffix, ".wrl");
/*********************************************************************************************************************************************************************************/

// atomic colors and radii, center of mass, setting user modified atomic scales and colors
atoms_fill_basic_settings(xyz, natoms, atoms);
if(centerofmass==1)
 {
 for(j=0;j<3;j++) centerofmasscoor[j] = 0.0; sumofmass = 0.0;
 for(i=0;i<natoms;i++) {for(j=0;j<3;j++) centerofmasscoor[j] += xyz[i].coor[j]*xyz[i].mass; sumofmass += xyz[i].mass;}
 for(j=0;j<3;j++) centerofmasscoor[j] /= sumofmass;
 for(i=0;i<natoms;i++) {for(j=0;j<3;j++) xyz[i].coor[j] -= centerofmasscoor[j];}
 }
atoms_change_atomscales(xyz, natoms, input);
atoms_change_atomcolors(xyz, natoms, input);

// connectivity matrix and its modifications specified by user
atoms_make_connectivitymatrix(xyz, natoms, input);

// reading gauges and making their settings
ngauges = 0; for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "addgauge"))!=NULL) ngauges++;
if(ngauges!=0)
 {
 if((gauge = (GAUGE *) malloc(ngauges*sizeof(GAUGE)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 atoms_make_gauge_settings(gauge, ngauges, xyz, natoms, nmolecules, bondradius, input);
 }

// checking common part of output except arrows (those are checked somewhere else, as well as meshes)
check_basic_settings(camspec, camera, scenerot, scenerotangle, camrot, camerarot, camrq,
                     camerarq, camsky, camerasky, camzoom, camerazoom, finish, finset);
check_xyz(xyz, natoms, nmolecules, hideatoms, atomcolor, bondcolor, atomradiusscale,
          lightatomradiusscale, bondradius, defatomcolor, defbondcolor);
check_gauges(xyz, natoms, gauge, ngauges);

// making arrows
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "arrow"))!=NULL) break; i++;
for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endarrow"))!=NULL) break;
if(input[i-1]!=NULL && (input[j]==NULL || j<=i)) {puts("\nProblems when reading arrows data...\n"); exit(1);}
if(input[i-1]!=NULL) narrows = j - i;
 else narrows = 0;
if(narrows!=0)
 {
 if((tyc = (TYCKA *) malloc(narrows*sizeof(TYCKA)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((kon = (KONUS *) malloc(narrows*sizeof(KONUS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 arrow_make_arrows(tyc, kon, narrows, xyz, natoms, nmolecules, convertxyzval, centerofmass, centerofmasscoor, arrowradius, coneradiusscale, input);
 }

// making vibrations
if(extfileformat==0 || extfileformat==2 || extfileformat==3) vibrations = 0;
if(extfileformat==1)
 {
 for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "vibration"))!=NULL) break; i++;
 for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endvibration"))!=NULL) break;
 if(input[i-1]!=NULL && (input[j]==NULL || j<=i)) {puts("\nProblems when reading vibrations data...\n"); exit(1);}
 if(input[i-1]!=NULL) vibrations = 1;
  else vibrations = 0;
 }
if(vibrations==1) vib_make_vibrations(xyz, natoms, nmolecules, centerofmass, centerofmasscoor, convertxyzval,
                                      hideatoms, atomcolor, bondcolor, atomradiusscale, lightatomradiusscale,
                                      bondradius, defatomcolor, defbondcolor, tyc, kon, narrows,
                                      gauge, ngauges, camspec, useangle, camera, scenerot, scenerotangle,
                                      camrot, camerarot, camrq, camerarq, camsky, camerasky, camzoom, camerazoom,
                                      handedness, finish, finset, input, extinput, extfileformat, outfilename,
                                      suffix, outformat, command, comend);

// making orbitals
if(extfileformat==0 || extfileformat==3) orbitals = 0;
if(extfileformat==1 || extfileformat==2)
 {
 for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "orbital"))!=NULL) break; i++;
 for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endorbital"))!=NULL) break;
 if(input[i-1]!=NULL && (input[j]==NULL || j<=i)) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
 if(input[i-1]!=NULL) orbitals = 1;
  else orbitals = 0;
 }
if(orbitals==1) orb_make_orbitals(xyz, natoms, nmolecules, centerofmass, centerofmasscoor, convertxyzval,
                                  hideatoms, atomcolor, bondcolor, atomradiusscale, lightatomradiusscale,
                                  bondradius, defatomcolor, defbondcolor, tyc, kon, narrows,
                                  gauge, ngauges, camspec, useangle, camera, scenerot, scenerotangle,
                                  camrot, camerarot, camrq, camerarq, camsky, camerasky, camzoom, camerazoom,
                                  handedness, finish, finset, input, extinput, extfileformat, outfilename,
                                  suffix, outformat, command, comend);

// making surfaces
if(extfileformat==0 || extfileformat==1 || extfileformat==2) surfaces = 0;
if(extfileformat==3)
 {
 for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "surface"))!=NULL) break; i++;
 for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endsurface"))!=NULL) break;
 if(input[i-1]!=NULL && (input[j]==NULL || j<=i)) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
 if(input[i-1]!=NULL) surfaces = 1;
  else surfaces = 0;
 }
if(surfaces==1) surf_make_surfaces(xyz, natoms, nmolecules, centerofmass, centerofmasscoor, convertxyzval,
                                   hideatoms, atomcolor, bondcolor, atomradiusscale, lightatomradiusscale,
                                   bondradius, defatomcolor, defbondcolor, tyc, kon, narrows,
                                   gauge, ngauges, camspec, useangle, camera, scenerot, scenerotangle,
                                   camrot, camerarot, camrq, camerarq, camsky, camerasky, camzoom, camerazoom,
                                   handedness, finish, finset, input, extinput, extfileformat, outfilename,
                                   suffix, outformat, command, comend);

// writing output file and possibly running povray
// this is done only if vibrations, orbitals and surfaces were not done
// (these functions produce their own outputs)
if(vibrations==0 && orbitals==0 && surfaces==0)
 {
 strcat(outfilename, suffix);
 if((outfile = fopen(outfilename, "w"))==NULL) {printf("\nOpening file %s failed!\n\n", outfilename); exit(1);}
 // common part of output
 outfile = write_basic_settings(xyz, natoms, tyc, kon, narrows, camspec, useangle, camera, convertxyzval, scenerot,
                                scenerotangle, camrot, camerarot, camrq, camerarq, camsky, camerasky, camzoom, camerazoom,
                                handedness, finish, finset, outfile, outformat);
 outfile = write_xyz(xyz, natoms, nmolecules, hideatoms, atomcolor, bondcolor, finish, atomradiusscale,
                     lightatomradiusscale, bondradius, defatomcolor, defbondcolor, outfile, outformat);
 outfile = write_gauges(xyz, natoms, gauge, ngauges, hideatoms, finish, atomradiusscale,
                        lightatomradiusscale, bondradius, outfile, outformat);
 outfile = write_arrows(tyc, kon, narrows, finish, outfile, outformat);
 if((fclose(outfile))==EOF) {printf("\nClosing file %s failed!\n\n", outfilename); exit(1);}
 // running povray
 if(outformat==1 && (read_firststrcmp(command, "norendering"))==NULL)
  {
  strcat(command, outfilename);
  strcat(command, comend);
  write_run_system_command(command);
  }
 }

// deallocate memory
if(extinput!=NULL) {for(i=0;extinput[i]!=NULL;i++) free(extinput[i]); free(extinput); extinput = NULL;}
for(i=0;input[i]!=NULL;i++) free(input[i]); free(input); input = NULL;
for(i=0;atoms[i]!=NULL;i++) free(atoms[i]); free(atoms); atoms = NULL;
for(i=0;i<natoms;i++) free(xyz[i].cone); free(xyz); xyz = NULL;
for(i=0;i<=nmolecules;i++) {free(defatomcolor[i]); free(defbondcolor[i]);}
free(defatomcolor); defatomcolor = NULL;
free(defbondcolor); defbondcolor = NULL;
free(tyc); tyc = NULL;
free(kon); kon = NULL;
free(ids); ids = NULL;
free(bondradius); bondradius = NULL;
free(arrowradius); arrowradius = NULL;
free(atomradiusscale); atomradiusscale = NULL;
free(lightatomradiusscale); lightatomradiusscale = NULL;
free(bondradiusscale); bondradiusscale = NULL;
free(arrowradiusscale); arrowradiusscale = NULL;
free(coneradiusscale); coneradiusscale = NULL;
free(gauge); gauge = NULL;

return 0;
}

