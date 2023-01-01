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
#include "charorb.h"
#include "charread.h"
#include "charcheck.h"
#include "charwrite.h"
#include "charcalc.h"
#include "charmesh.h"

/*********************************************************************************************************************************************************************************/
// functions for orbitals
// rename Cs symmetry labels because povray doesn't support ' and " in filenames anymore
void orb_label_renamecs(char *outlabel, const char * const inplabel)
{
int i;
static int renamecs = 0;

strcpy(outlabel, inplabel);

i = 0; while(outlabel[i] != '\0') i++; i--;

if(outlabel[i] == '\'')
 {
 outlabel[i] = '\0';
 strcat(outlabel, "p");
 renamecs++;
 }
else if(outlabel[i] == '\"')
 {
 outlabel[i] = '\0';
 strcat(outlabel, "pp");
 renamecs++;
 }

if(renamecs == 1) puts("Replacing \' by p and \" by pp in orbital-symmetry labels for charmol output filenames because of povray's limited functionality...");
}

// cutting orbital information into 4 strings (num1, sym1, num2, sym2)
char *orb_cut_strings(char *final, char *orig)
{
int i, j, c, o1n, o2n;
char o1s[20], o2s[20], str1[50], str2[50];

if((strstr(orig, "-"))==NULL)
 {
 if((c = sscanf(orig, "%d%s", &o1n, o1s))!=2)
  {
  if(c==1) strcpy(o1s, "order");
   else {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  }
 o2n = o1n;
 strcpy(o2s, o1s);
 }
 else
  {
  for(i=0,j=0; orig[i]!='-' ;i++,j++) str1[j] = orig[i]; str1[j] = '\0';
  for(i++,j=0; orig[i]!='\0' ;i++,j++) str2[j] = orig[i]; str2[j] = '\0';
  if((c = sscanf(str1, "%d%s", &o1n, o1s))==2)
   {
   if((c = sscanf(str2, "%d%s", &o2n, o2s))!=2)
    {
    if(c==1) {puts("\nSpecify symmetry for either both of orbital labels or none of orbital ordering numbers in orbital range...\n"); exit(1);}
     else {puts("\nProblems when reading orbitals data...\n"); exit(1);}
    }
   }
   else
    {
    if(c==1)
     {
     if((c = sscanf(str2, "%d%s", &o2n, o2s))>1)
      {puts("\nSpecify symmetry for either both of orbital labels or none of orbital ordering numbers in orbital range...\n"); exit(1);}
     if(c<1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
     strcpy(o1s, "order"); strcpy(o2s, o1s);
     }
     else {puts("\nProblems when reading orbitals data...\n"); exit(1);}
    }
  }

if((strcmp(o1s, o2s))!=0) {puts("\nOrbital ranges can be specified only for orbital ordering numbers or orbital labels within one symmetry...\n"); exit(1);}
if(o1n<o2n) sprintf(final, "%d  %d  %s", o1n, o2n, o1s);
 else sprintf(final, "%d  %d  %s", o2n, o1n, o1s);

return final;
}

// process options from a line
int orb_process_options(char *line, ORBITAL *orb, MOS *mo, int calconly, int norbs, int global, ORBCOMB *comb, int combine)
{
int i, j, k, l, c, tookmore, o1n, o2n, nactorbs, na, nb, naall, nball, ncomb, ival, gencon, screen, longname, style;
int isovalspec, screenfactorspec, rgbposspec, rgbnegspec, transparencyspec, genconspec, screenspec, gridsizescalespec, griddensityspec, longnamespec, stylespec;
static int combind = -1;
char *string, *token, os[20], *help, keyhelp[50], griddensity[20];
double isoval, screenfactor, rgbpos[3], rgbneg[3], transparency, gridsizescale;
typedef struct {int order; char label[20]; int spin; int combind;} ACTORBINFO;
ACTORBINFO *actorb;

const char keywords[] = " alpha beta isovalue generalcontraction prescreening prescreeningfactor gridsizescale griddensity positivecolor negativecolor transparency longfilenames style combineorbitals ";

help = NULL;
if((help = (char *) malloc(100*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
actorb = NULL;
if((actorb = (ACTORBINFO *) malloc(1*sizeof(ACTORBINFO)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
actorb[0].spin = 0;

na = 0; nb = 0; naall = 0; nball = 0; ncomb = 0;
gencon = 0; screen = 0; longname = 0; style = 0;
isovalspec = 0; screenfactorspec = 0; rgbposspec = 0; rgbnegspec = 0; transparencyspec = 0;
genconspec = 0; screenspec = 0; gridsizescalespec = 0; griddensityspec = 0; longnamespec = 0; stylespec = 0;

string = line;
token = read_next_token(string);
string = NULL;

while(token!=NULL)
 {
 strcpy(keyhelp, " "); strcat(keyhelp, token); strcat(keyhelp, " ");
 if((strstr(keywords, keyhelp))==NULL) {printf("\nUnknown keyword '%s' in orbital section...\n\n", token); exit(1);}
 tookmore = 0;
 if((strcmp(token, "isovalue"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((c = sscanf(token, "%lf", &isoval))!=1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  isovalspec = 1;
  }
 if((strcmp(token, "generalcontraction"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((strcmp(token, "on"))!=0 && (strcmp(token, "off"))!=0) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((strcmp(token, "off"))==0) gencon = 0;
   else gencon = 1;
  genconspec = 1;
  }
 if((strcmp(token, "prescreening"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((strcmp(token, "on"))!=0 && (strcmp(token, "off"))!=0) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((strcmp(token, "off"))==0) screen = 0;
   else screen = 1;
  screenspec = 1;
  }
 if((strcmp(token, "prescreeningfactor"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((c = sscanf(token, "%lf", &screenfactor))!=1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  screenfactorspec = 1;
  }
 if((strcmp(token, "gridsizescale"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((c = sscanf(token, "%lf", &gridsizescale))!=1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  gridsizescalespec = 1;
  }
 if((strcmp(token, "griddensity"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((strcmp(token, "verylow"))!=0 && (strcmp(token, "low"))!=0 && (strcmp(token, "medium"))!=0 &&
     (strcmp(token, "high"))!=0 && (strcmp(token, "veryhigh"))!=0) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  strcpy(griddensity, token);
  griddensityspec = 1;
  }
 if((strcmp(token, "positivecolor"))==0)
  {
  for(i=0;i<3;i++)
   {
   token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
   if((c = sscanf(token, "%lf", &rgbpos[i]))!=1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
   }
  rgbposspec = 1;
  }
 if((strcmp(token, "negativecolor"))==0)
  {
  for(i=0;i<3;i++)
   {
   token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
   if((c = sscanf(token, "%lf", &rgbneg[i]))!=1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
   }
  rgbnegspec = 1;
  }
 if((strcmp(token, "transparency"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((c = sscanf(token, "%lf", &transparency))!=1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  transparencyspec = 1;
  }
 if((strcmp(token, "longfilenames"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((strcmp(token, "on"))!=0 && (strcmp(token, "off"))!=0) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((strcmp(token, "on"))==0) longname = 1;
   else longname = 0;
  longnamespec = 1;
  }
 if((strcmp(token, "style"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  style = 0;
  if((strcmp(token, "solid"))==0)      style = 1;
  if((strcmp(token, "wireframe"))==0)  style = 2;
  if(style==0) {puts("\nProblems when reading surfaces data...\n"); exit(1);}
  stylespec = 1;
  }
 if((strcmp(token, "combineorbitals"))==0)
  {
  ncomb++;
  if(calconly==0) combind++;
  }
 if(token!=NULL && (strcmp(token, "alpha"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((strcmp(token, "all"))==0)
   {
   tookmore = 0;
   if(naall!=1)
    {
    na = 0; for(i=0;mo[i].coef!=NULL;i++) if(mo[i].spin==1) na++; naall = 1;
    if(na==0) {puts("\nNo alpha orbitals found in external file...\n"); exit(1);}
    if(calconly==0)
     {
     i = 0; while(actorb[i].spin!=0) i++;
     if((actorb = (ACTORBINFO *) realloc(actorb, (na+i+1)*sizeof(ACTORBINFO)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
     for(j=0,k=0;j<na;j++)
      {
      while(mo[k].spin!=1) k++;
      strcpy(actorb[j+i].label, mo[k].label); k++;
      actorb[j+i].spin = 1;
      actorb[j+i].combind = combind;
      }
     actorb[j+i].spin = 0;
     }
    }
   }
   else
    {
    tookmore = 1;
    while(token!=NULL && (c = sscanf(token, "%d", &ival))==1)
     {
     if(naall!=1)
      {
      help = orb_cut_strings(help, token);
      if((c = sscanf(help, "%d%d%s", &o1n, &o2n, os))!=3) {puts("\nProblems in orb_cut_strings function...\n"); exit(1);}
      na += o2n - o1n + 1;
      if(calconly==0)
       {
       i = 0; while(actorb[i].spin!=0) i++;
       if((actorb = (ACTORBINFO *) realloc(actorb, ((o2n-o1n+1)+i+1)*sizeof(ACTORBINFO)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
       for(j=0;(j+o1n)<=o2n;j++)
        {
        sprintf(actorb[j+i].label, "%d", j+o1n);
        strcat(actorb[j+i].label, os);
        actorb[j+i].spin = 1;
        actorb[j+i].combind = combind;
        }
       actorb[j+i].spin = 0;
       }
      }
     token = read_next_token(string);
     }
    }
  }
 if(token!=NULL && (strcmp(token, "beta"))==0)
  {
  token = read_next_token(string); if(token==NULL) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
  if((strcmp(token, "all"))==0)
   {
   tookmore = 0;
   if(nball!=1)
    {
    nb = 0; for(i=0;mo[i].coef!=NULL;i++) if(mo[i].spin==-1) nb++; nball = 1;
    if(nb==0) {puts("\nNo beta orbitals found in external file...\n"); exit(1);}
    if(calconly==0)
     {
     i = 0; while(actorb[i].spin!=0) i++;
     if((actorb = (ACTORBINFO *) realloc(actorb, (nb+i+1)*sizeof(ACTORBINFO)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
     for(j=0,k=0;j<nb;j++)
      {
      while(mo[k].spin!=-1) k++;
      strcpy(actorb[j+i].label, mo[k].label); k++;
      actorb[j+i].spin = -1;
      actorb[j+i].combind = combind;
      }
     actorb[j+i].spin = 0;
     }
    }
   }
   else
    {
    tookmore = 1;
    while(token!=NULL && (c = sscanf(token, "%d", &ival))==1)
     {
     if(nball!=1)
      {
      help = orb_cut_strings(help, token);
      if((c = sscanf(help, "%d%d%s", &o1n, &o2n, os))!=3) {puts("\nProblems in orb_cut_strings function...\n"); exit(1);}
      nb += o2n - o1n + 1;
      if(calconly==0)
       {
       i = 0; while(actorb[i].spin!=0) i++;
       if((actorb = (ACTORBINFO *) realloc(actorb, ((o2n-o1n+1)+i+1)*sizeof(ACTORBINFO)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
       for(j=0;(j+o1n)<=o2n;j++)
        {
        sprintf(actorb[j+i].label, "%d", j+o1n);
        strcat(actorb[j+i].label, os);
        actorb[j+i].spin = -1;
        actorb[j+i].combind = combind;
        }
       actorb[j+i].spin = 0;
       }
      }
     token = read_next_token(string);
     }
    }
  }
 if(token!=NULL && tookmore!=1) token = read_next_token(string);
 }

if(combine==0)
 {
 nactorbs = na + nb;
 if(nactorbs==0 && global==0) {puts("\nNo orbitals to calculate were specified on line containing keyword alpha or beta...\n"); exit(1);}
 if(calconly==0)
  {
  if(global==1)
   {
   for(i=0;i<norbs;i++)
    {
    if(isovalspec==1) orb[i].isoval = isoval;
    if(genconspec==1) orb[i].gencon = gencon;
    if(screenspec==1) orb[i].screen = screen;
    if(screenfactorspec==1) orb[i].screenfactor = screenfactor;
    if(gridsizescalespec==1) orb[i].gridsizescale = gridsizescale;
    if(griddensityspec==1) strcpy(orb[i].griddensity, griddensity);
    if(rgbposspec==1) for(k=0;k<3;k++) orb[i].rgbpos[k] = rgbpos[k];
    if(rgbnegspec==1) for(k=0;k<3;k++) orb[i].rgbneg[k] = rgbneg[k];
    if(transparencyspec==1) orb[i].transparency = transparency;
    if(longnamespec==1) orb[i].longname = longname;
    if(stylespec==1) orb[i].style = style;
    }
   }
   else
    {
    i = 0; while(orb[i].spin!=0) i++;
    for(j=0;actorb[j].spin!=0;j++,i++)
     {
     orb[i].spin = actorb[j].spin;
     if((strstr(actorb[j].label, "order"))==NULL)
      {
      strcpy(orb[i].label, actorb[j].label);
      k = 0; l = 0; while(mo[k].coef!=NULL) {if(mo[k].spin==orb[i].spin) {l++; if((strcmp(mo[k].label, orb[i].label))==0) break;} k++;}
      if(mo[k].coef==NULL)
       {
       if(orb[i].spin==1) {printf("\nOrbital %s alpha was not found in external file...\n\n", orb[i].label); exit(1);}
        else {printf("\nOrbital %s beta was not found in external file...\n\n", orb[i].label); exit(1);}
       }
      orb[i].order = l;
      }
      else
       {
       if((c = sscanf(actorb[j].label, "%d", &orb[i].order))!=1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
       k = 0; l = 0; while(mo[k].coef!=NULL) {if(mo[k].spin==orb[i].spin) {l++; if(orb[i].order==l) break;} k++;}
       if(mo[k].coef==NULL)
        {
        if(orb[i].spin==1) {printf("\nOrbital %d alpha was not found in external file...\n\n", orb[i].order); exit(1);}
         else {printf("\nOrbital %d beta was not found in external file...\n\n", orb[i].order); exit(1);}
        }
       strcpy(orb[i].label, mo[k].label);
       }
     k = 0; while(mo[k].coef!=NULL) {if(orb[i].spin==mo[k].spin && (strcmp(orb[i].label, mo[k].label))==0) break; else k++;}
     if(mo[k].coef==NULL)
      {
      if(orb[i].spin==1) {printf("\nOrbital %s alpha was not found in external file...\n\n", orb[i].label); exit(1);}
       else {printf("\nOrbital %s beta was not found in external file...\n\n", orb[i].label); exit(1);}
      }
      else orb[i].coef = mo[k].coef;
     if(isovalspec==1) orb[i].isoval = isoval;
     if(genconspec==1) orb[i].gencon = gencon;
     if(screenspec==1) orb[i].screen = screen;
     if(screenfactorspec==1) orb[i].screenfactor = screenfactor;
     if(gridsizescalespec==1) orb[i].gridsizescale = gridsizescale;
     if(griddensityspec==1) strcpy(orb[i].griddensity, griddensity);
     if(rgbposspec==1) for(k=0;k<3;k++) orb[i].rgbpos[k] = rgbpos[k];
     if(rgbnegspec==1) for(k=0;k<3;k++) orb[i].rgbneg[k] = rgbneg[k];
     if(transparencyspec==1) orb[i].transparency = transparency;
     if(longnamespec==1) orb[i].longname = longname;
     if(stylespec==1) orb[i].style = style;
     }
    }
  }
 }
 else
  {
  if(calconly==0)
   {
   nactorbs = na + nb;
   if(nactorbs==0) {puts("\nNo orbitals to combine were specified...\n"); exit(1);}
   for(i=0;actorb[i].spin!=0;i++)
    {
    j = actorb[i].combind;
    if(j==-1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
    comb[j].norb++;
    if((comb[j].orbs = (int *) realloc(comb[j].orbs, comb[j].norb*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
    if((strstr(actorb[i].label, "order"))==NULL)
     {
     k = 0; while(k<norbs) {if(orb[k].spin==actorb[i].spin && (strcmp(orb[k].label, actorb[i].label))==0) break; else k++;}
     }
     else
      {
      if((c = sscanf(actorb[i].label, "%d", &actorb[i].order))!=1) {puts("\nProblems when reading orbitals data...\n"); exit(1);}
      k = 0; while(k<norbs) {if(orb[k].spin==actorb[i].spin && orb[k].order==actorb[i].order) break; else k++;}
      }
    if(k==norbs) {puts("\nOrbitals to combine has to be first specified by separate alpha or beta keyword...\n"); exit(1);}
     else comb[j].orbs[comb[j].norb - 1] = k;
    }   
   }
   else nactorbs = ncomb;
  }

free(actorb); actorb = NULL;
free(help); help = NULL;

return nactorbs;
}

/*********************************************************************************************************************************************************************************/
// main function for orbitals
// making orbitals
void orb_make_orbitals(ATOMS *xyz, int natoms, int nmolecules, int centerofmass, double *centerofmasscoor, double convertxyzval,
                       int hideatoms, int atomcolor, int bondcolor, double *atomradiusscale, double *lightatomradiusscale,
                       double *bondradius, double **defatomcolor, double **defbondcolor, TYCKA *tyc, KONUS *kon, int narrows,
                       GAUGE *gauge, int ngauges, int camspec, int useangle, double *camera, int scenerot, double scenerotangle,
                       int camrot, double *camerarot, int camrq, double *camerarq, int camsky, double *camerasky, int camzoom,
                       double camerazoom, int handedness, int finish, double *finset, char **input, char **extinput,
                       int extfileformat, char *outfilename, char *suffix, int outformat, char *command, char *comend)
{
int i, j, k, l, o, nshell, nbasis, norbs, ncombs, maxcombs, gpoints[3], totpoints, maxdegen, maxnprim, gencon;
char actfilename[500], actcommand[550], namehelp[300], labelhelp[100];
double *gridx, *gridy, *gridz, *gridval, *gridexp, gmin[3], val;
FILE *outfile;
BASIS *gto;
MOS *mo;
ORBITAL *orb;
ORBCOMB *comb;
MESH mesh;

// initialization
gto = NULL; mo = NULL; orb = NULL; comb = NULL; gridx = NULL; gridy = NULL; gridz = NULL; gridval = NULL; gridexp = NULL;

// reading GTO basis set and MO coefficients
// specification of format of external files: 1 - molden, 2 - fchk
nshell = 0; nbasis = 0;
switch(extfileformat)
 {
 case 1:
  gto = read_gto_molden(xyz, natoms, extinput);
  while(gto[nshell].atom!=-1) {nbasis += gto[nshell].degen; nshell++;}
  mo = read_mo_molden(nbasis, extinput);
  break;
 case 2:
  gto = read_gto_fchk(xyz, natoms, extinput);
  while(gto[nshell].atom!=-1) {nbasis += gto[nshell].degen; nshell++;}
  mo = read_mo_fchk(nbasis, extinput);
  break;
 }

// determining number of orbitals to calculate and preparing for settings
for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "orbital"))!=NULL) break; i++;
for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endorbital"))!=NULL) break;
norbs = 0;
for(k=i;k<j;k++)
 {
 if(((strstr(input[k], "alpha"))!=NULL || (strstr(input[k], "beta"))!=NULL)
    && (strstr(input[k], "combineorbitals"))==NULL) norbs += orb_process_options(input[k], orb, mo, 1, norbs, 0, comb, 0);
 }
if(norbs==0) {puts("\nNo orbital specification found in 'orbital' section...\n"); exit(1);}
if((orb = (ORBITAL *) malloc(norbs*sizeof(ORBITAL)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}

// determining number of orbital combinations
ncombs = 0;
for(k=i;k<j;k++)
 {
 if((strstr(input[k], "combineorbitals"))!=NULL) ncombs += orb_process_options(input[k], orb, mo, 1, norbs, 0, comb, 1);
 }
if((comb = (ORBCOMB *) malloc(ncombs*sizeof(ORBCOMB)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}

/*********************************************************************************************************************************************************************************/
// VARIOUS SETTINGS IMPORTANT TO USERS

// default settings for orbitals
for(k=0;k<norbs;k++)
 {
 orb[k].isoval = 0.05;
 orb[k].gencon = 1;
 orb[k].screen = 1;
 orb[k].screenfactor = 1E-3;
 orb[k].gridsizescale = 1.0;
 strcpy(orb[k].griddensity, "medium");
 orb[k].transparency = 0.0;
 orb[k].rgbpos[0] = 0.0; orb[k].rgbpos[1] = 0.0; orb[k].rgbpos[2] = 1.0;
 orb[k].rgbneg[0] = 1.0; orb[k].rgbneg[1] = 0.0; orb[k].rgbneg[2] = 0.0;
 orb[k].coef = NULL;
 orb[k].order = 0;
 orb[k].spin = 0;
 orb[k].longname = 0;
 orb[k].style = 1;
 }

// modification of global settings for orbitals specified by user
for(k=i;k<j;k++)
 {
 if((strstr(input[k], "alpha"))==NULL && (strstr(input[k], "beta"))==NULL
    && (strstr(input[k], "combineorbitals"))==NULL) orb_process_options(input[k], orb, mo, 0, norbs, 1, comb, 0);
 }

// modification of settings for particular orbitals specified by user
// orbital labels, spin and coef pointers are filled during this
for(k=i;k<j;k++)
 {
 if(((strstr(input[k], "alpha"))!=NULL || (strstr(input[k], "beta"))!=NULL)
    && (strstr(input[k], "combineorbitals"))==NULL) orb_process_options(input[k], orb, mo, 0, norbs, 0, comb, 0);
 }

// finishing the settings for orbitals
for(k=0;k<norbs;k++)
 {
 orb[k].screenval = fabs(orb[k].screenfactor*orb[k].isoval);
 if((strcmp(orb[k].griddensity, "verylow"))==0)  for(l=0;l<3;l++) orb[k].griddist[l] = 0.60;
 if((strcmp(orb[k].griddensity, "low"))==0)      for(l=0;l<3;l++) orb[k].griddist[l] = 0.35;
 if((strcmp(orb[k].griddensity, "medium"))==0)   for(l=0;l<3;l++) orb[k].griddist[l] = 0.18;
 if((strcmp(orb[k].griddensity, "high"))==0)     for(l=0;l<3;l++) orb[k].griddist[l] = 0.10;
 if((strcmp(orb[k].griddensity, "veryhigh"))==0) for(l=0;l<3;l++) orb[k].griddist[l] = 0.05;
 }

// default settings for combinations of orbitals
for(k=0;k<ncombs;k++)
 {
 comb[k].file = NULL;
 comb[k].norb = 0;
 comb[k].orbs = NULL;
 }

// modification of settings for particular combinations of orbitals specified by user
for(k=i;k<j;k++)
 {
 if((strstr(input[k], "combineorbitals"))!=NULL) orb_process_options(input[k], orb, mo, 0, norbs, 0, comb, 1);
 }
/*********************************************************************************************************************************************************************************/

// checking number of orbital combinations according to maximum number of open files FOPEN_MAX from stdio.h
#ifdef FOPEN_MAX
 maxcombs = (int) FOPEN_MAX/2;
#else
 maxcombs = 8;
#endif
if(ncombs>maxcombs) {printf("\nNumber of orbital combinations is limited to %d...\n\n", maxcombs); exit(1);}

// checking settings of orbitals for hard problems
for(o=0;o<norbs;o++)
 {
 check_mesh(orb[o].rgbpos, orb[o].transparency, orb[o].isoval, 0);
 check_mesh(orb[o].rgbneg, orb[o].transparency, -orb[o].isoval, 0);
 }

// determining the highest degeneration and number of primitives of used shells
maxdegen = 0; maxnprim = 0;
for(i=0;i<nshell;i++)
 {
 if(gto[i].degen>maxdegen) maxdegen = gto[i].degen;
 if(gto[i].nprim>maxnprim) maxnprim = gto[i].nprim;
 }

// checking whether it is generally-contracted basis set
gencon = calc_check_gencon(natoms, nshell, maxdegen, gto);
if(gencon==0)
 {
 for(i=0;i<norbs;i++) if(orb[i].gencon==1) orb[i].gencon = 0;
 }

// preparing orbital combinations
for(i=0;i<ncombs;i++)
 {
 strcpy(comb[i].filename, outfilename);
 strcat(comb[i].filename, ".comb");
 for(j=0;j<comb[i].norb;j++) if(orb[comb[i].orbs[j]].spin==1) {strcat(comb[i].filename, ".alpha"); break;}
 while(j<comb[i].norb)
  {
  if(orb[comb[i].orbs[j]].spin==1)
   {
   sprintf(namehelp, ".%d.", orb[comb[i].orbs[j]].order);
   strcat(comb[i].filename, namehelp);
   //strcat(comb[i].filename, orb[comb[i].orbs[j]].label); // Cs symmetry label hack because of povray
   orb_label_renamecs(labelhelp, orb[comb[i].orbs[j]].label);
   strcat(comb[i].filename, labelhelp);
   }
  j++;
  }
 for(j=0;j<comb[i].norb;j++) if(orb[comb[i].orbs[j]].spin==-1) {strcat(comb[i].filename, ".beta"); break;}
 while(j<comb[i].norb)
  {
  if(orb[comb[i].orbs[j]].spin==-1)
   {
   sprintf(namehelp, ".%d.", orb[comb[i].orbs[j]].order);
   strcat(comb[i].filename, namehelp);
   //strcat(comb[i].filename, orb[comb[i].orbs[j]].label); // Cs symmetry label hack because of povray
   orb_label_renamecs(labelhelp, orb[comb[i].orbs[j]].label);
   strcat(comb[i].filename, labelhelp);
   }
  j++;
  }
 strcat(comb[i].filename, suffix);
 if((comb[i].file = fopen(comb[i].filename, "w"))==NULL) {printf("\nOpening file %s failed!\n\n", comb[i].filename); exit(1);}
 // common part of output
 comb[i].file = write_basic_settings(xyz, natoms, tyc, kon, narrows, camspec, useangle, camera, convertxyzval, scenerot,
                                     scenerotangle, camrot, camerarot, camrq, camerarq, camsky, camerasky, camzoom, camerazoom,
                                     handedness, finish, finset, comb[i].file, outformat);
 comb[i].file = write_xyz(xyz, natoms, nmolecules, hideatoms, atomcolor, bondcolor, finish, atomradiusscale,
                          lightatomradiusscale, bondradius, defatomcolor, defbondcolor, comb[i].file, outformat);
 comb[i].file = write_gauges(xyz, natoms, gauge, ngauges, hideatoms, finish, atomradiusscale,
                             lightatomradiusscale, bondradius, comb[i].file, outformat);
 comb[i].file = write_arrows(tyc, kon, narrows, finish, comb[i].file, outformat);
 }

// calculation and depiction of all required orbitals in a loop
// specification of format of external files: 1 - molden, 2 - fchk
for(o=0;o<norbs;o++)
 {
 // determining the grid for calculation and its generation
 for(i=0;i<3;i++)
  {
  gmin[i] = xyz[0].coor[i];
  for(j=1;j<natoms;j++) if(xyz[j].coor[i]<gmin[i]) gmin[i] = xyz[j].coor[i];
  val = xyz[0].coor[i];
  for(j=1;j<natoms;j++) if(xyz[j].coor[i]>val) val = xyz[j].coor[i];
  gmin[i] -= 5.0;
  val += 5.0;
  gpoints[i] = (int) (val - gmin[i])*orb[o].gridsizescale/orb[o].griddist[i];
  if((gpoints[i]%2)==0) gpoints[i] += 1;
  gmin[i] += (val - gmin[i])/2.0 - orb[o].griddist[i]*(gpoints[i] - 1)/2.0;
  }
 totpoints = gpoints[0]*gpoints[1]*gpoints[2];
 if((gridx = (double *) realloc(gridx, totpoints*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((gridy = (double *) realloc(gridy, totpoints*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((gridz = (double *) realloc(gridz, totpoints*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((gridval = (double *) realloc(gridval, totpoints*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
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
    gridx[l] = gmin[0] + i*orb[o].griddist[0];
    gridy[l] = gmin[1] + j*orb[o].griddist[1];
    gridz[l] = gmin[2] + k*orb[o].griddist[2];
    }
   }
  }
 // checking whether we can afford genereally-contracted scheme due to higher memory requirements (if yes, memory for gridexp is reallocated)
 if(orb[o].gencon==1)
  {
  if((gridexp = (double *) realloc(gridexp, (maxnprim*totpoints)*sizeof(double)))==NULL)
   {
   orb[o].gencon = 0; puts("Switching off using general contraction because of too high memory requirements...");
   }
  }
 // predeterming screening radii (squared)
 if(orb[o].screen==1) calc_screenr2_values(gto, natoms, nshell, maxdegen, orb[o].screenval, orb[o].gencon);
 // calculating values in grid points with or without using general contraction
 switch(extfileformat)
  {
  case 1:
   if(orb[o].gencon==1) calc_grid_gencon_molden(gto, natoms, nshell, maxdegen, totpoints, orb[o].screen, gridx, gridy, gridz, gridval, orb[o].coef, gridexp, orb[o].gencon);
    else calc_grid_segcon_molden(gto, natoms, nshell, maxdegen, totpoints, orb[o].screen, gridx, gridy, gridz, gridval, orb[o].coef, gridexp, orb[o].gencon);
   break;
  case 2:
   if(orb[o].gencon==1) calc_grid_gencon_fchk(gto, natoms, nshell, maxdegen, totpoints, orb[o].screen, gridx, gridy, gridz, gridval, orb[o].coef, gridexp, orb[o].gencon);
    else calc_grid_segcon_fchk(gto, natoms, nshell, maxdegen, totpoints, orb[o].screen, gridx, gridy, gridz, gridval, orb[o].coef, gridexp, orb[o].gencon);
   break;
  }
 // writing output
 strcpy(actfilename, outfilename);
 if(orb[o].spin==1) strcat(actfilename, ".alpha.");
  else strcat(actfilename, ".beta.");
 sprintf(namehelp, "%d.", orb[o].order);
 strcat(actfilename, namehelp);
 //strcat(actfilename, orb[o].label); // Cs symmetry label hack because of povray
 orb_label_renamecs(labelhelp, orb[o].label);
 strcat(actfilename, labelhelp);
 if(orb[o].longname==1)
  {
  if(orb[o].screen==1) sprintf(namehelp, ".iso-%f.gridscale-%.3f.griddens-%s.screenfac-%f.posrgb-%.3f-%.3f-%.3f.negrgb-%.3f-%.3f-%.3f.transp-%.3f",
                               orb[o].isoval, orb[o].gridsizescale, orb[o].griddensity, orb[o].screenfactor, orb[o].rgbpos[0], orb[o].rgbpos[1],
                               orb[o].rgbpos[2], orb[o].rgbneg[0], orb[o].rgbneg[1], orb[o].rgbneg[2], orb[o].transparency);
   else sprintf(namehelp, ".iso-%f.gridscale-%.3f.griddens-%s.posrgb-%.3f-%.3f-%.3f.negrgb-%.3f-%.3f-%.3f.transp-%.3f",
                orb[o].isoval, orb[o].gridsizescale, orb[o].griddensity, orb[o].rgbpos[0], orb[o].rgbpos[1],
                orb[o].rgbpos[2], orb[o].rgbneg[0], orb[o].rgbneg[1], orb[o].rgbneg[2], orb[o].transparency);
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
 mesh = mesh_generate_mesh(gpoints, orb[o].griddist, gridx, gridy, gridz, gridval, NULL, orb[o].isoval);
 switch(orb[o].style)
  {
  case 1: outfile = write_mesh_solid(mesh, orb[o].rgbpos, orb[o].transparency, orb[o].isoval, finish, outfile, outformat); break;
  case 2: outfile = write_mesh_wireframe(mesh, orb[o].rgbpos, orb[o].transparency, orb[o].isoval, finish, outfile, outformat); break;
  }
 for(i=0;i<ncombs;i++)
  {
  for(j=0;j<comb[i].norb;j++)
   {
   if(comb[i].orbs[j]==o)
    {
    switch(orb[o].style)
     {
     case 1: comb[i].file = write_mesh_solid(mesh, orb[o].rgbpos, orb[o].transparency, orb[o].isoval, finish, comb[i].file, outformat); break;
     case 2: comb[i].file = write_mesh_wireframe(mesh, orb[o].rgbpos, orb[o].transparency, orb[o].isoval, finish, comb[i].file, outformat); break;
     }
    }
   }
  }
 mesh = mesh_clean_mesh(mesh);
 mesh = mesh_generate_mesh(gpoints, orb[o].griddist, gridx, gridy, gridz, gridval, NULL, -orb[o].isoval);
 switch(orb[o].style)
  {
  case 1: outfile = write_mesh_solid(mesh, orb[o].rgbneg, orb[o].transparency, -orb[o].isoval, finish, outfile, outformat); break;
  case 2: outfile = write_mesh_wireframe(mesh, orb[o].rgbneg, orb[o].transparency, -orb[o].isoval, finish, outfile, outformat); break;
  }
 for(i=0;i<ncombs;i++)
  {
  for(j=0;j<comb[i].norb;j++)
   {
   if(comb[i].orbs[j]==o)
    {
    switch(orb[o].style)
     {
     case 1: comb[i].file = write_mesh_solid(mesh, orb[o].rgbneg, orb[o].transparency, -orb[o].isoval, finish, comb[i].file, outformat); break;
     case 2: comb[i].file = write_mesh_wireframe(mesh, orb[o].rgbneg, orb[o].transparency, -orb[o].isoval, finish, comb[i].file, outformat); break;
     }
    }
   }
  }
 mesh = mesh_clean_mesh(mesh);
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

// finishing orbital combinations
for(i=0;i<ncombs;i++)
 {
 if((fclose(comb[i].file))==EOF) {printf("\nClosing file %s failed!\n\n", comb[i].filename); exit(1);}
 // running povray
 if(outformat==1 && (read_firststrcmp(command, "norendering"))==NULL)
  {
  strcpy(actcommand, command);
  strcat(actcommand, comb[i].filename);
  strcat(actcommand, comend);
  write_run_system_command(actcommand);
  }
 }

// deallocate memory
for(i=0;i<nshell;i++) {free(gto[i].expon); free(gto[i].contr); if(gto[i].degen==4) free(gto[i].contr2);} free(gto); gto = NULL;
for(i=0;mo[i].coef!=NULL;i++) free(mo[i].coef); free(mo); mo = NULL;
for(i=0;i<ncombs;i++) free(comb[i].orbs); free(comb); comb = NULL;
free(orb); orb = NULL;
free(gridx); gridx = NULL;
free(gridy); gridy = NULL;
free(gridz); gridz = NULL;
free(gridval); gridval = NULL;
free(gridexp); gridexp = NULL;
}

