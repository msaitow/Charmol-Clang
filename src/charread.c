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
#include "charread.h"

/*********************************************************************************************************************************************************************************/
// functions for reading files to char matrices
// reading input file to input char matrix
char **read_file_to_char(char *filename)
{
FILE *f;
int rad, sl, nsl, c;
char **input;

input = NULL;
rad=0; sl=0;
if((f = fopen(filename, "r"))==NULL) {printf("\nOpening file %s failed!\n\n", filename); exit(1);}
if((c=fgetc(f))==EOF) {printf("\nFile %s is empty!\n\n", filename); fclose(f); exit(1);}
do
 {
 nsl = 100;
 if((input = (char **) realloc(input, (rad+1)*sizeof(char *)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
 if((input[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
 if(c=='\n') input[rad][0]='\0';
  else
   {
   input[rad][0]=(char)c;
   for(sl=1;(c=fgetc(f))!='\n' && c!=EOF;sl++)
    {
    input[rad][sl]=(char)c;
    if(sl==(nsl-3))
     {
     nsl += 50;
     if((input[rad] = (char *) realloc(input[rad], nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
     }
    }
   input[rad][sl]='\0';
   }
 rad++;
 }
while((c=fgetc(f))!=EOF);
if((input = (char **) realloc(input, (rad+1)*sizeof(char *)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
input[rad] = NULL;
fclose(f);

input = read_repair_fortran_single_precision(input, filename);

return input;
}

// reading xyz part of cube file to input char matrix
char **read_xyzpartofcubefile_to_char(char *filename)
{
FILE *f;
int rad, sl, nsl, c, d, natoms;
char **input;

input = NULL;
rad=0; sl=0; natoms=0;
if((f = fopen(filename, "r"))==NULL) {printf("\nOpening file %s failed!\n\n", filename); exit(1);}
if((c=fgetc(f))==EOF) {printf("\nFile %s is empty!\n\n", filename); fclose(f); exit(1);}
do
 {
 nsl = 100;
 if((input = (char **) realloc(input, (rad+1)*sizeof(char *)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
 if((input[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
 if(c=='\n') input[rad][0]='\0';
  else
   {
   input[rad][0]=(char)c;
   for(sl=1;(c=fgetc(f))!='\n' && c!=EOF;sl++)
    {
    input[rad][sl]=(char)c;
    if(sl==(nsl-3))
     {
     nsl += 50;
     if((input[rad] = (char *) realloc(input[rad], nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
     }
    }
   input[rad][sl]='\0';
   if(rad==2)
    {
    if((d = sscanf(input[rad], "%d", &natoms))!=1) {puts("\nProblems when reading xyz part of cube file...\n"); exit(1);}
    if(natoms<0) natoms = -natoms;
    natoms += 15;
    }
   }
 rad++;
 if(rad==natoms) break;
 }
while((c=fgetc(f))!=EOF);
if((input = (char **) realloc(input, (rad+1)*sizeof(char *)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
input[rad] = NULL;
fclose(f);

input = read_repair_fortran_single_precision(input, filename);

return input;
}

// reading input file to input char matrix without blank lines and comment lines (beginning with #)
char **read_file_to_char_withoutcomments(char *filename)
{
FILE *f;
int rad, sl, nsl, c;
char **input;

input = NULL;
rad=0; sl=0;
if((f = fopen(filename, "r"))==NULL) {printf("\nOpening file %s failed!\n\n", filename); exit(1);}
if((c=fgetc(f))==EOF) {printf("\nFile %s is empty!\n\n", filename); fclose(f); exit(1);}
do
 {
 if(c=='#' || c=='\n')
  {
  if(c=='#') {
    while((c=fgetc(f))!='\n' && c!=EOF)
      if(c==EOF) break;
  }
  }
  else
   {
   nsl = 100;
   if((input = (char **) realloc(input, (rad+1)*sizeof(char *)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
   if((input[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
   if(c=='\n') input[rad][0]='\0';
    else
     {
     input[rad][0]=(char)c;
     for(sl=1;(c=fgetc(f))!='\n' && c!=EOF;sl++)
      {
      input[rad][sl]=(char)c;
      if(sl==(nsl-3))
       {
       nsl += 50;
       if((input[rad] = (char *) realloc(input[rad], nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
       }
      }
     input[rad][sl]='\0';
     }
   rad++;
   }
 }
while((c=fgetc(f))!=EOF);
if((input = (char **) realloc(input, (rad+1)*sizeof(char *)))==NULL) {puts("\nMemory allocation failed!\n"); fclose(f); exit(1);}
input[rad] = NULL;
fclose(f);

input = read_repair_fortran_single_precision(input, filename);

return input;
}

/*********************************************************************************************************************************************************************************/
// functions for modifications of the char matrices (input data)
// repairing fortran single precision exponents (replacing 'D' by 'E')
char **read_repair_fortran_single_precision(char **input, char *filename)
{
int i, j, k, repaired, number, ndig, nsig, ndot;
char *line;

repaired = 0;

for(i=0;input[i]!=NULL;i++)
 {
 line = input[i];
 if(strstr(line, "file")==NULL) // rather leave lines with file specifications untouched
  {
  for(j=1;line[j]!='\0';j++) // number can not begin with 'd' ('e')
   {
   if(line[j]=='d' || line[j]=='D') // now check if it is part of a number
    {
    number = 1;
    // check to the right
    ndig = 0; nsig = 0;
    for(k=j+1;line[k]!='\0' && line[k]!=' ' && line[k]!='\t' && line[k]!=',';k++)
     {
     if(line[k]>='0' && line[k]<='9') ndig++;
      else if(line[k]=='+' || line[k]=='-') nsig++;
       else {number = 0; break;}
     }
    if(number==1)
     {
     if(ndig<1 || nsig>1) number = 0;
      else if(nsig==1 && line[j+1]!='+' && line[j+1]!='-') number = 0;
     }
    // check to the left
    if(number==1)
     {
     ndig = 0; nsig = 0; ndot = 0;
     for(k=j-1;k>=0 && line[k]!=' ' && line[k]!='\t' && line[k]!=',';k--)
      {
      if(line[k]>='0' && line[k]<='9') ndig++;
       else if(line[k]=='+' || line[k]=='-') nsig++;
        else if(line[k]=='.') ndot++;
         else {number = 0; break;}
      }
     if(number==1)
      {
      if(ndig<1 || nsig>1 || ndot>1) number = 0;
       else if(nsig==1 && line[k+1]!='+' && line[k+1]!='-') number = 0;
      }
     }
    // if it is a number, replace 'D' by 'E'
    if(number==1)
     {
     line[j] = 'E';
     repaired = 1;
     }
    }
   }
  }
 }

if(repaired==1) printf("Fortran single precision exponents from '%s' file replaced with C format...\n", filename);

return input;
}

/*********************************************************************************************************************************************************************************/
// functions helping with interpreting strings
// taking next token in a string
char *read_next_token(char *string)
{
int i;
char *tok;
static char *str = NULL, chr = '\0';

if(str==NULL && string==NULL) return NULL;
if(string!=NULL) str = string;
 else
  {
  if(chr!='\0') {str[0] = chr; chr = '\0';}
  }

i = 0;
while((str[i]==' ' || str[i]=='\t') && str[i]!='\0') i++;
if(str[i]=='\0') return NULL;
 else
  {
  tok = &str[i];
  i = 0;
  while(tok[i]!=' ' && tok[i]!='\t' && tok[i]!='\0') i++;
  if(tok[i]!='\0') {chr = tok[i]; tok[i] = '\0';}
  str = &tok[i];
  }

return tok;
}

// compare first tokens in two strings
char *read_firststrcmp(char *str, char *tok)
{
int i, k;
char *isit;

if(str==NULL || tok==NULL) return NULL;

i = 0; k = 0;
while((str[i]==' ' || str[i]=='\t') && str[i]!='\0') i++;
while((tok[k]==' ' || tok[k]=='\t') && tok[k]!='\0') k++;
if(str[i]=='\0' || tok[k]=='\0') return NULL;

while(tok[k]==str[i] && tok[k]!=' ' && tok[k]!='\t' && tok[k]!='\0') {i++; k++;}

if((tok[k]==' ' || tok[k]=='\t' || tok[k]=='\0') && (str[i]==' ' || str[i]=='\t' || str[i]=='\0')) isit = &tok[k];
 else isit = NULL;

return isit;
}

/*********************************************************************************************************************************************************************************/
// functions for reading XYZ coordinates from char matrices
// reading xyz from charmol input
ATOMS *read_xyz_input(char **input)
{
int i, j, k, c, natoms, nmolecules;
ATOMS *tmpxyz;

tmpxyz = NULL;

for(i=0;input[i]!=NULL;i++) if((read_firststrcmp(input[i], "xyz"))!=NULL) break; i++;
for(j=0;input[j]!=NULL;j++) if((read_firststrcmp(input[j], "endxyz"))!=NULL) break;
if(input[i-1]==NULL || input[j]==NULL || j<=i) {puts("\nGeometry not found in charmol input file...\n"); exit(1);}
 else natoms = j - i;

if((tmpxyz = (ATOMS *) malloc((natoms+1)*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
k = 0;
if((c = sscanf(input[k+i], "%s%lf%lf%lf%d", tmpxyz[k].symb, &tmpxyz[k].coor[0], &tmpxyz[k].coor[1], &tmpxyz[k].coor[2], &tmpxyz[k].molid))==5) nmolecules = 1;
 else {nmolecules = 0; tmpxyz[k].molid = 0; if(c!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}}
for(k=1;k<natoms;k++)
 {
 if(nmolecules==0)
  {
  if((c = sscanf(input[k+i], "%s%lf%lf%lf", tmpxyz[k].symb, &tmpxyz[k].coor[0], &tmpxyz[k].coor[1], &tmpxyz[k].coor[2]))!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}
  tmpxyz[k].molid = 0;
  }
  else
   {
   if((c = sscanf(input[k+i], "%s%lf%lf%lf%d", tmpxyz[k].symb, &tmpxyz[k].coor[0], &tmpxyz[k].coor[1], &tmpxyz[k].coor[2], &tmpxyz[k].molid))!=5)
    {puts("\nProblems when reading xyz data...\nIf you have specified more molecules, make sure that you have also specified mid for each atom.\n"); exit(1);}
   }
 }
tmpxyz[k].atomnum = -1;

return tmpxyz;
}

// reading xyz from molden
ATOMS *read_xyz_molden(char **input)
{
int i, j, c, natoms, nmolecules, format;
ATOMS *tmpxyz;

tmpxyz = NULL;

format = 0;
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "[ATOMS]"))!=NULL || (strstr(input[i], "[Atoms]"))!=NULL || (strstr(input[i], "[atoms]"))!=NULL) break;
if(input[i]!=NULL) format = 1;
 else
  {
  if((c = sscanf(input[0], "%d", &natoms))==1) format = 2;
   else
    {
    for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "[FR-COORD]"))!=NULL || (strstr(input[i], "[Fr-coord]"))!=NULL || (strstr(input[i], "[fr-coord]"))!=NULL) break;
    if(input[i]!=NULL) format = 3;
    }
  }
if(format==0) {puts("\nGeometry not found in molden file...\n"); exit(1);}

switch(format)
 {
 case 1: // ATOMS
  for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "[ATOMS]"))!=NULL || (strstr(input[i], "[Atoms]"))!=NULL || (strstr(input[i], "[atoms]"))!=NULL) break; i++;
  j = 0;
  if((tmpxyz = (ATOMS *) malloc(2*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
  if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
  if((c = sscanf(input[i], "%*s%*d%d%lf%lf%lf%d", &tmpxyz[j].atomnum, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2], &tmpxyz[j].molid))==5) nmolecules = 1;
   else {nmolecules = 0; tmpxyz[j].molid = 0; if(c!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}}
  tmpxyz[j].symb[0] = '\0';
  j++; i++;
  if(input[i]!=NULL)
   {
   if(nmolecules==0)
    {
    while((c = sscanf(input[i], "%*s%*d%d%lf%lf%lf", &tmpxyz[j].atomnum, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2]))==4)
     {
     tmpxyz[j].symb[0] = '\0';
     tmpxyz[j].molid = 0;
     j++; i++;
     if((tmpxyz = (ATOMS *) realloc(tmpxyz, (j+1)*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
     if(input[i]==NULL) break;
     }
    }
    else
     {
     while((c = sscanf(input[i], "%*s%*d%d%lf%lf%lf%d", &tmpxyz[j].atomnum, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2], &tmpxyz[j].molid))==5)
      {
      tmpxyz[j].symb[0] = '\0';
      j++; i++;
      if((tmpxyz = (ATOMS *) realloc(tmpxyz, (j+1)*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
      if(input[i]==NULL) break;
      }
     if(c==4) puts("Make sure that you have specified mid for each atom. It is possible, that some atoms were omitted during read of xyz from molden file...");
     }
   }
  break;
 case 2: // simple
  if(natoms<1) {puts("\nNumber of atoms smaller than 1 does not make sense...\n"); exit(1);}
  if((tmpxyz = (ATOMS *) malloc((natoms+1)*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
  i = 2; j = 0;
  if(input[i]==NULL || input[i-1]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
  if((c = sscanf(input[i], "%s%lf%lf%lf%d", tmpxyz[j].symb, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2], &tmpxyz[j].molid))==5) nmolecules = 1;
   else {nmolecules = 0; tmpxyz[j].molid = 0; if(c!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}}
  for(j++,i++;j<natoms;j++,i++)
   {
   if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
   if(nmolecules==0)
    {
    if((c = sscanf(input[i], "%s%lf%lf%lf", tmpxyz[j].symb, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2]))!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}
    tmpxyz[j].molid = 0;
    }
    else
     {
     if((c = sscanf(input[i], "%s%lf%lf%lf%d", tmpxyz[j].symb, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2], &tmpxyz[j].molid))!=5)
      {puts("\nProblems when reading xyz data...\nIf you have specified more molecules, make sure that you have also specified mid for each atom.\n"); exit(1);}
     }
   }
  break;
 case 3: // FR-COORD
  for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "[FR-COORD]"))!=NULL || (strstr(input[i], "[Fr-coord]"))!=NULL || (strstr(input[i], "[fr-coord]"))!=NULL) break; i++;
  j = 0;
  if((tmpxyz = (ATOMS *) malloc(2*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
  if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
  if((c = sscanf(input[i], "%s%lf%lf%lf%d", tmpxyz[j].symb, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2], &tmpxyz[j].molid))==5) nmolecules = 1;
   else {nmolecules = 0; tmpxyz[j].molid = 0; if(c!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}}
  j++; i++;
  if(input[i]!=NULL)
   {
   if(nmolecules==0)
    {
    while((c = sscanf(input[i], "%s%lf%lf%lf", tmpxyz[j].symb, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2]))==4)
     {
     tmpxyz[j].molid = 0;
     j++; i++;
     if((tmpxyz = (ATOMS *) realloc(tmpxyz, (j+1)*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
     if(input[i]==NULL) break;
     }
    }
    else
     {
     while((c = sscanf(input[i], "%s%lf%lf%lf%d", tmpxyz[j].symb, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2], &tmpxyz[j].molid))==5)
      {
      j++; i++;
      if((tmpxyz = (ATOMS *) realloc(tmpxyz, (j+1)*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
      if(input[i]==NULL) break;
      }
     if(c==4) puts("Make sure that you have specified mid for each atom. It is possible, that some atoms were omitted during read of xyz from molden file...");
     }
   }
  break;
 }
tmpxyz[j].atomnum = -1;

return tmpxyz;
}

// reading xyz from fchk
ATOMS *read_xyz_fchk(char **input)
{
int i, j, k, l, c, natoms, nmolecules, ncoor;
char *string, *token;
ATOMS *tmpxyz;

tmpxyz = NULL;

for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Number of atoms"))!=NULL) break;
if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%d", &natoms))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
if(natoms<1) {puts("\nNumber of atoms smaller than 1 does not make sense...\n"); exit(1);}
if((tmpxyz = (ATOMS *) malloc((natoms+1)*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Atomic numbers"))!=NULL) break;
if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%d", &j))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
if(j!=natoms) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
k = 0;
while(k<natoms)
 {
 i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%d", &tmpxyz[k].atomnum))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
 k++;
 while((token = read_next_token(string))!=NULL)
  {
  if(k<natoms)
   {
   if((c = sscanf(token, "%d", &tmpxyz[k].atomnum))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
   k++;
   }
  }
 }
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Molecular ids"))!=NULL) break;
if(input[i]==NULL) nmolecules = 0;
 else nmolecules = 1;
if(nmolecules==1)
 {
 if((c = sscanf(input[i], "%*s%*s%*s%*s%d", &j))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
 if(j!=natoms) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 k = 0;
 while(k<natoms)
  {
  i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
  string = input[i];
  if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
  string = NULL;
  if((c = sscanf(token, "%d", &tmpxyz[k].molid))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
  k++;
  while((token = read_next_token(string))!=NULL)
   {
   if(k<natoms)
    {
    if((c = sscanf(token, "%d", &tmpxyz[k].molid))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
    k++;
    }
   }
  }
 }
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Current cartesian coordinates"))!=NULL || (strstr(input[i], "Current Cartesian coordinates"))!=NULL) break;
if(input[i]==NULL) {puts("\nGeometry not found in fchk file...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%*s%d", &ncoor))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
if(ncoor!=3*natoms) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
j = 0; k = 0; l = 0;
while(j<ncoor)
 {
 i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%lf", &tmpxyz[k].coor[l]))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
 j++; l++; if(l==3) {l = 0; k++;}
 while((token = read_next_token(string))!=NULL)
  {
  if(j<ncoor)
   {
   if((c = sscanf(token, "%lf", &tmpxyz[k].coor[l]))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
   j++; l++; if(l==3) {l = 0; k++;}
   }
  }
 }
tmpxyz[k].atomnum = -1;
for(i=0;tmpxyz[i].atomnum!=-1;i++) tmpxyz[i].symb[0] = '\0';
if(nmolecules==0) {for(i=0;tmpxyz[i].atomnum!=-1;i++) tmpxyz[i].molid = 0;}

return tmpxyz;
}

// reading xyz from cube
ATOMS *read_xyz_cube(char **input)
{
int i, j, c, natoms, nmolecules;
ATOMS *tmpxyz;

tmpxyz = NULL;

if(input[0]==NULL || input[1]==NULL || input[2]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
i = 2;
if((c = sscanf(input[i], "%d", &natoms))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
if(natoms<0) natoms = -natoms; // negative number of atoms for orbital cube file
if(natoms<1) {puts("\nNumber of atoms smaller than 1 does not make sense...\n"); exit(1);}
if((tmpxyz = (ATOMS *) malloc((natoms+1)*sizeof(ATOMS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
if(input[3]==NULL || input[4]==NULL || input[5]==NULL || input[6]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
i = 6; j = 0;
if((c = sscanf(input[i], "%d%*f%lf%lf%lf%d", &tmpxyz[j].atomnum, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2], &tmpxyz[j].molid))==5) nmolecules = 1;
 else {nmolecules = 0; tmpxyz[j].molid = 0; if(c!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}}
for(j++,i++;j<natoms;j++,i++)
 {
 if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 if(nmolecules==0)
  {
  if((c = sscanf(input[i], "%d%*f%lf%lf%lf", &tmpxyz[j].atomnum, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2]))!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}
  tmpxyz[j].molid = 0;
  }
  else
   {
   if((c = sscanf(input[i], "%d%*f%lf%lf%lf%d", &tmpxyz[j].atomnum, &tmpxyz[j].coor[0], &tmpxyz[j].coor[1], &tmpxyz[j].coor[2], &tmpxyz[j].molid))!=5)
    {puts("\nProblems when reading xyz data...\nIf you have specified more molecules, make sure that you have also specified mid for each atom.\n"); exit(1);}
   }
 }
tmpxyz[j].atomnum = -1;
for(i=0;i<natoms;i++) tmpxyz[i].symb[0] = '\0';

return tmpxyz;
}

/*********************************************************************************************************************************************************************************/
// functions for reading GTOs from char matrices
// reading GTOs from molden (cartesian d, f and g functions are used in molden by default, GTOs in atomic units)
// and preparing constants (as much of the constant as possible together with contraction coefficient is put to contr or contr2)
BASIS *read_gto_molden(ATOMS *xyz, int natoms, char **input)
{
int nd, nf, ng, i, j, k, l, c, test;
BASIS *tmpgto;

tmpgto = NULL;

nd = 6; nf = 10; ng = 15;
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "5D"))!=NULL) {nd = 5; break;}
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "7F"))!=NULL) {nf = 7; break;}
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "9G"))!=NULL) {ng = 9; break;}
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "[GTO]"))!=NULL || (strstr(input[i], "[Gto]"))!=NULL || (strstr(input[i], "[gto]"))!=NULL) break; i++;
if(input[i-1]==NULL) {puts("\nNo GTO basis found in molden file...\n"); exit(1);}
if((tmpgto = (BASIS *) malloc(1*sizeof(BASIS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
for(j=0,k=0;j<natoms;j++)
 {
 if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
 if((c = sscanf(input[i], "%d", &tmpgto[k].atom))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);} i++;
 if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
 while((c = sscanf(input[i], "%s%d", tmpgto[k].type, &tmpgto[k].nprim))==2)
  {
  i++;
  if((tmpgto[k].expon = (double *) malloc(tmpgto[k].nprim*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
  if((tmpgto[k].contr = (double *) malloc(tmpgto[k].nprim*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
  if((strstr(tmpgto[k].type, "sp"))!=NULL)
   {
   if((tmpgto[k].contr2 = (double *) malloc(tmpgto[k].nprim*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
   for(l=0;l<tmpgto[k].nprim;l++)
    {
    if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
    if((c = sscanf(input[i], "%lf%lf%lf", &tmpgto[k].expon[l], &tmpgto[k].contr[l], &tmpgto[k].contr2[l]))!=3)
     {puts("\nProblems when reading GTO data...\n"); exit(1);} i++;
    }
   }
   else
    {
    for(l=0;l<tmpgto[k].nprim;l++)
     {
     if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
     if((c = sscanf(input[i], "%lf%lf", &tmpgto[k].expon[l], &tmpgto[k].contr[l]))!=2) {puts("\nProblems when reading GTO data...\n"); exit(1);} i++;
     }
    }
  for(l=0;l<3;l++) tmpgto[k].coor[l] = xyz[tmpgto[k].atom-1].coor[l];
  k++;
  if((tmpgto = (BASIS *) realloc(tmpgto, (k+1)*sizeof(BASIS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
  tmpgto[k].atom = tmpgto[k-1].atom;
  if(input[i]==NULL) break;
  }
 // there does not have to be always empty line (but should be there)
 if((c = sscanf(input[i], "%d", &test))!=1) i++;
 if((input[i-1]==NULL || input[i]==NULL) && j<(natoms-1)) {puts("\nGTO basis was not read for all atoms. Molden file seems to be corrupted...\n"); exit(1);}
 }
tmpgto[k].atom = -1;
for(i=0;tmpgto[i].atom!=-1;i++)
 {
 tmpgto[i].degen = 0;
 // orbital s (whole constant)
 if((strstr(tmpgto[i].type, "s"))!=NULL && (strstr(tmpgto[i].type, "sp"))==NULL)
  {
  tmpgto[i].degen = 1;
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75);
  }
 // orbital p (whole constant)
 if((strstr(tmpgto[i].type, "p"))!=NULL && (strstr(tmpgto[i].type, "sp"))==NULL)
  {
  tmpgto[i].degen = 3;
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*2.0*sqrt(tmpgto[i].expon[j]);
  }
 // orbital sp (whole constant)
 if((strstr(tmpgto[i].type, "sp"))!=NULL)
  {
  tmpgto[i].degen = 4;
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75);
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr2[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*2.0*sqrt(tmpgto[i].expon[j]);
  }
 // orbital d (partial constant: (2*alpha/pi)^(3/4)*4*alpha)
 if((strstr(tmpgto[i].type, "d"))!=NULL)
  {
  tmpgto[i].degen = nd;
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*4.0*tmpgto[i].expon[j];
  }
 // orbital f (partial constant: (2*alpha/pi)^(3/4)*(4*alpha)^(3/2))
 if((strstr(tmpgto[i].type, "f"))!=NULL)
  {
  tmpgto[i].degen = nf;
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*pow(4.0*tmpgto[i].expon[j], 1.5);
  }
 // orbital g (partial constant: (2*alpha/pi)^(3/4)*(4*alpha)^2)
 if((strstr(tmpgto[i].type, "g"))!=NULL)
  {
  tmpgto[i].degen = ng;
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*pow(4.0*tmpgto[i].expon[j], 2.0);
  }
 if(tmpgto[i].degen == 0) {puts("\nMolden format supports angular momenta only up to 'g'...\n"); exit(1);}
 }

return tmpgto;
}

// reading GTOs from fchk (GTOs in atomic units)
// and preparing constants (as much of the constant as possible together with contraction coefficient is put to contr or contr2)
BASIS *read_gto_fchk(ATOMS *xyz, int natoms, char **input)
{
int i, j, k, l, c, nbasis, nshell, nexpcon, spis;
char *string, *token;
BASIS *tmpgto;

tmpgto = NULL;

for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Number of basis functions"))!=NULL) break;
if(input[i]==NULL) {puts("\nNo GTO basis found in fchk file...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%*s%d", &nbasis))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Number of contracted shells"))!=NULL) break;
if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%*s%d", &nshell))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
if((tmpgto = (BASIS *) malloc((nshell+1)*sizeof(BASIS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Shell types"))!=NULL) break;
if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%d", &j))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
if(j!=nshell) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
k = 0;
while(k<nshell)
 {
 i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%d", &tmpgto[k].degen))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
 k++;
 while((token = read_next_token(string))!=NULL)
  {
  if(k<nshell)
   {
   if((c = sscanf(token, "%d", &tmpgto[k].degen))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
   k++;
   }
  }
 }
tmpgto[k].atom = -1;
for(k=0;tmpgto[k].atom!=-1;k++)
 {
 if(abs(tmpgto[k].degen)>5) {puts("\nSorry, charmol program can not handle higher than 'h' angular momentum yet...\n"); exit(1);}
 switch(tmpgto[k].degen)
  {
  // Shell types in fchk: 0=s, 1=p, -1=sp, 2=6d, -2=5d, 3=10f, -3=7f, 4=15g, -4=9g, 5=21h, -5=11h
  case  0: strcpy(tmpgto[k].type, "s");  tmpgto[k].degen =  1; break;
  case  1: strcpy(tmpgto[k].type, "p");  tmpgto[k].degen =  3; break;
  case -1: strcpy(tmpgto[k].type, "sp"); tmpgto[k].degen =  4; break;
  case  2: strcpy(tmpgto[k].type, "d");  tmpgto[k].degen =  6; break;
  case -2: strcpy(tmpgto[k].type, "d");  tmpgto[k].degen =  5; break;
  case  3: strcpy(tmpgto[k].type, "f");  tmpgto[k].degen = 10; break;
  case -3: strcpy(tmpgto[k].type, "f");  tmpgto[k].degen =  7; break;
  case  4: strcpy(tmpgto[k].type, "g");  tmpgto[k].degen = 15; break;
  case -4: strcpy(tmpgto[k].type, "g");  tmpgto[k].degen =  9; break;
  case  5: strcpy(tmpgto[k].type, "h");  tmpgto[k].degen = 21; break;
  case -5: strcpy(tmpgto[k].type, "h");  tmpgto[k].degen = 11; break;
  }
 }
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Number of primitives per shell"))!=NULL) break;
if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%*s%*s%*s%d", &j))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
if(j!=nshell) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
k = 0;
while(k<nshell)
 {
 i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%d", &tmpgto[k].nprim))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
 k++;
 while((token = read_next_token(string))!=NULL)
  {
  if(k<nshell)
   {
   if((c = sscanf(token, "%d", &tmpgto[k].nprim))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
   k++;
   }
  }
 }
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Shell to atom map"))!=NULL) break;
if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%*s%*s%d", &j))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
if(j!=nshell) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
k = 0;
while(k<nshell)
 {
 i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%d", &tmpgto[k].atom))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
 k++;
 while((token = read_next_token(string))!=NULL)
  {
  if(k<nshell)
   {
   if((c = sscanf(token, "%d", &tmpgto[k].atom))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
   k++;
   }
  }
 }
nexpcon = 0; spis = 0;
for(k=0;tmpgto[k].atom!=-1;k++)
 {
 for(l=0;l<3;l++) tmpgto[k].coor[l] = xyz[tmpgto[k].atom-1].coor[l];
 nexpcon += tmpgto[k].nprim;
 if((tmpgto[k].expon = (double *) malloc(tmpgto[k].nprim*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((tmpgto[k].contr = (double *) malloc(tmpgto[k].nprim*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if(tmpgto[k].degen==4)
  {
  spis = 1;
  if((tmpgto[k].contr2 = (double *) malloc(tmpgto[k].nprim*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
  }
 }
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Primitive exponents"))!=NULL) break;
if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%d", &j))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
if(j!=nexpcon) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
j = 0; k = 0; l = 0;
while(j<nexpcon)
 {
 i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%lf", &tmpgto[k].expon[l]))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
 j++; l++; if(l==tmpgto[k].nprim) {l = 0; k++;}
 while((token = read_next_token(string))!=NULL)
  {
  if(j<nexpcon)
   {
   if((c = sscanf(token, "%lf", &tmpgto[k].expon[l]))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
   j++; l++; if(l==tmpgto[k].nprim) {l = 0; k++;}
   }
  }
 }
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Contraction coefficients"))!=NULL && (strstr(input[i], "P(S=P)"))==NULL) break;
if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%d", &j))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
if(j!=nexpcon) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
j = 0; k = 0; l = 0;
while(j<nexpcon)
 {
 i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%lf", &tmpgto[k].contr[l]))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
 j++; l++; if(l==tmpgto[k].nprim) {l = 0; k++;}
 while((token = read_next_token(string))!=NULL)
  {
  if(j<nexpcon)
   {
   if((c = sscanf(token, "%lf", &tmpgto[k].contr[l]))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
   j++; l++; if(l==tmpgto[k].nprim) {l = 0; k++;}
   }
  }
 }
if(spis==1)
 {
 for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "P(S=P) Contraction coefficients"))!=NULL) break;
 if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 if((c = sscanf(input[i], "%*s%*s%*s%*s%*s%d", &j))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}
 if(j!=nexpcon) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 j = 0; k = 0; l = 0;
 while(j<nexpcon)
  {
  i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
  string = input[i];
  if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
  string = NULL;
  if(tmpgto[k].degen==4) {if((c = sscanf(token, "%lf", &tmpgto[k].contr2[l]))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}}
  j++; l++; if(l==tmpgto[k].nprim) {l = 0; k++;}
  while((token = read_next_token(string))!=NULL)
   {
   if(j<nexpcon)
    {
    if(tmpgto[k].degen==4) {if((c = sscanf(token, "%lf", &tmpgto[k].contr2[l]))!=1) {puts("\nProblems when reading GTO data...\n"); exit(1);}}
    j++; l++; if(l==tmpgto[k].nprim) {l = 0; k++;}
    }
   }
  }
 }
for(i=0;tmpgto[i].atom!=-1;i++)
 {
 // orbital s (whole constant)
 if((strstr(tmpgto[i].type, "s"))!=NULL && (strstr(tmpgto[i].type, "sp"))==NULL)
  {
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75);
  }
 // orbital p (whole constant)
 if((strstr(tmpgto[i].type, "p"))!=NULL && (strstr(tmpgto[i].type, "sp"))==NULL)
  {
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*2.0*sqrt(tmpgto[i].expon[j]);
  }
 // orbital sp (whole constant)
 if((strstr(tmpgto[i].type, "sp"))!=NULL)
  {
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75);
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr2[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*2.0*sqrt(tmpgto[i].expon[j]);
  }
 // orbital d (partial constant: (2*alpha/pi)^(3/4)*4*alpha)
 if((strstr(tmpgto[i].type, "d"))!=NULL)
  {
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*4.0*tmpgto[i].expon[j];
  }
 // orbital f (partial constant: (2*alpha/pi)^(3/4)*(4*alpha)^(3/2))
 if((strstr(tmpgto[i].type, "f"))!=NULL)
  {
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*pow(4.0*tmpgto[i].expon[j], 1.5);
  }
 // orbital g (partial constant: (2*alpha/pi)^(3/4)*(4*alpha)^2)
 if((strstr(tmpgto[i].type, "g"))!=NULL)
  {
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*pow(4.0*tmpgto[i].expon[j], 2.0);
  }
 // orbital h (partial constant: (2*alpha/pi)^(3/4)*(4*alpha)^(5/2))
 if((strstr(tmpgto[i].type, "h"))!=NULL)
  {
  for(j=0;j<tmpgto[i].nprim;j++) tmpgto[i].contr[j] *= pow(2.0*tmpgto[i].expon[j]/M_PI, 0.75)*pow(4.0*tmpgto[i].expon[j], 2.5);
  }
 }

return tmpgto;
}

/*********************************************************************************************************************************************************************************/
// functions for reading MOs from char matrices
// reading MOs from molden
MOS *read_mo_molden(int nbasis, char **input)
{
int i, j, l, k, c, alpha, beta, repaired;
char *str;
MOS *tmpmo;

tmpmo = NULL;

for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "[MO]"))!=NULL || (strstr(input[i], "[Mo]"))!=NULL || (strstr(input[i], "[mo]"))!=NULL) break; i++;
if(input[i-1]==NULL) {puts("\nNo MOs found in molden file...\n"); exit(1);}
if((tmpmo = (MOS *) malloc((2*nbasis+1)*sizeof(MOS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
j = 0; alpha = 0; beta = 0;
while(input[i]!=NULL && j<(2*nbasis) && ((strstr(input[i], "Sym"))!=NULL || (strstr(input[i], "Spin"))!=NULL ||
                                         (strstr(input[i], "Ene"))!=NULL || (strstr(input[i], "Occup"))!=NULL))
 {
 if((tmpmo[j].coef = (double *) malloc(nbasis*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 if((strstr(input[i], "Sym"))==NULL && (strstr(input[i], "Spin"))==NULL &&
    (strstr(input[i], "Ene"))==NULL && (strstr(input[i], "Occup"))==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
 tmpmo[j].label[0] = '\0';
 tmpmo[j].spin = 0;
 for(l=0;l<4;l++,i++)
  {
  if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
  if((strstr(input[i], "Sym"))==NULL && (strstr(input[i], "Spin"))==NULL &&
     (strstr(input[i], "Ene"))==NULL && (strstr(input[i], "Occup"))==NULL) break;
  if((strstr(input[i], "Sym"))!=NULL)
   {
   if((c = sscanf(input[i], "%*s%s", tmpmo[j].label))!=1)
    {
    k = 0; while(input[i][k]!='\0' && input[i][k]!='=') k++;
    if(input[i][k]=='\0') {puts("\nProblems when reading MOs data...\n"); exit(1);}
     else
      {
      k++; str = &input[i][k];
      if((c = sscanf(str, "%s", tmpmo[j].label))!=1) {puts("\nProblems when reading MOs data...\n"); exit(1);}
      }
    }
   }
  if((strstr(input[i], "Spin"))!=NULL)
   {
   if((strstr(input[i], "Alpha"))!=NULL) {tmpmo[j].spin = 1; alpha++;}
    else {tmpmo[j].spin = -1; beta++;}
   }
  }
 if(tmpmo[j].spin==0) {puts("\nMolden file seems to be corrupted...\nSpin of MOs has to be specified!\n"); exit(1);}
 if(tmpmo[j].label[0]=='\0')
  {
  if(tmpmo[j].spin==1) sprintf(tmpmo[j].label, "%da", alpha);
   else sprintf(tmpmo[j].label, "%da", beta);
  }
 for(l=0;l<nbasis;l++,i++)
  {
  if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
  if((c = sscanf(input[i], "%*d%lf", &tmpmo[j].coef[l]))!=1) {puts("\nProblems when reading MOs data...\n"); exit(1);}
  }
 j++;
 }
tmpmo[j].coef = NULL;
for(i=0;tmpmo[i].coef!=NULL;i++)
 {
 j = 0; while(tmpmo[i].label[j]>='0' && tmpmo[i].label[j]<='9' && tmpmo[i].label[j]!='\0') j++;
 if(tmpmo[i].label[j]=='\0') {tmpmo[i].label[j] = 'a'; tmpmo[i].label[j+1] = '\0';}
 }

repaired = 0;
alpha = 1; beta = 1;
for(i=0;tmpmo[i].coef!=NULL;i++)
 {
 if(tmpmo[i].spin==1) {if(alpha>1 && tmpmo[i].label[0]=='1' && tmpmo[i].label[1]=='a' && tmpmo[i].label[2]=='\0') {sprintf(tmpmo[i].label, "%da", alpha); repaired = 1;} alpha++;}
  else                {if( beta>1 && tmpmo[i].label[0]=='1' && tmpmo[i].label[1]=='a' && tmpmo[i].label[2]=='\0') {sprintf(tmpmo[i].label, "%da",  beta); repaired = 1;}  beta++;}
 }
if(repaired==1) printf("ORCA's molden molecular orbital labels were repaired...\n");

return tmpmo;
}

// reading MOs from fchk
MOS *read_mo_fchk(int nbasis, char **input)
{
int i, j, k, l, c, order, nalpha, nalnum, nbeta, nbenum;
char *string, *token;
MOS *tmpmo;

tmpmo = NULL;

for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Alpha MO coefficients"))!=NULL) break;
if(input[i]==NULL) {puts("\nMOs not found in fchk file...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%*s%*s%*s%*s%d", &nalnum))!=1) {puts("\nProblems when reading MOs data...\n"); exit(1);}
if((nalnum%nbasis)!=0) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
nalpha = nalnum/nbasis;
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Beta MO coefficients"))!=NULL) break;
if(input[i]==NULL) {nbeta = 0; nbenum = 0;}
 else
  {
  if((c = sscanf(input[i], "%*s%*s%*s%*s%*s%d", &nbenum))!=1) {puts("\nProblems when reading MOs data...\n"); exit(1);}
  if((nbenum%nbasis)!=0) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
  nbeta = nbenum/nbasis;
  }
if((tmpmo = (MOS *) malloc((nalpha+nbeta+1)*sizeof(MOS)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
for(i=0;i<(nalpha+nbeta);i++) if((tmpmo[i].coef = (double *) malloc(nbasis*sizeof(double)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Alpha MO coefficients"))!=NULL) break;
j = 0; k = 0; l = 0; order = 1;
while(j<nalnum)
 {
 i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%lf", &tmpmo[k].coef[l]))!=1) {puts("\nProblems when reading MOs data...\n"); exit(1);}
 j++; l++; if(l==nbasis) {tmpmo[k].spin = 1; sprintf(tmpmo[k].label, "%d", order); l = 0; k++; order++;}
 while((token = read_next_token(string))!=NULL)
  {
  if(j<nalnum)
   {
   if((c = sscanf(token, "%lf", &tmpmo[k].coef[l]))!=1) {puts("\nProblems when reading MOs data...\n"); exit(1);}
   j++; l++; if(l==nbasis) {tmpmo[k].spin = 1; sprintf(tmpmo[k].label, "%d", order); l = 0; k++; order++;}
   }
  }
 }
if(nbenum!=0)
 {
 for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "Beta MO coefficients"))!=NULL) break;
 j = 0; order = 1;
 while(j<nbenum)
  {
  i++; if(input[i]==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
  string = input[i];
  if((token = read_next_token(string))==NULL) {puts("\nFchk file seems to be corrupted...\n"); exit(1);}
  string = NULL;
  if((c = sscanf(token, "%lf", &tmpmo[k].coef[l]))!=1) {puts("\nProblems when reading MOs data...\n"); exit(1);}
  j++; l++; if(l==nbasis) {tmpmo[k].spin = -1; sprintf(tmpmo[k].label, "%d", order); l = 0; k++; order++;}
  while((token = read_next_token(string))!=NULL)
   {
   if(j<nbenum)
    {
    if((c = sscanf(token, "%lf", &tmpmo[k].coef[l]))!=1) {puts("\nProblems when reading MOs data...\n"); exit(1);}
    j++; l++; if(l==nbasis) {tmpmo[k].spin = -1; sprintf(tmpmo[k].label, "%d", order); l = 0; k++; order++;}
    }
   }
  }
 }
tmpmo[k].coef = NULL;
for(i=0;tmpmo[i].coef!=NULL;i++)
 {
 j = 0; while(tmpmo[i].label[j]>='0' && tmpmo[i].label[j]<='9' && tmpmo[i].label[j]!='\0') j++;
 if(tmpmo[i].label[j]=='\0') {tmpmo[i].label[j] = 'a'; tmpmo[i].label[j+1] = '\0';}
 }

return tmpmo;
}

/*********************************************************************************************************************************************************************************/
// functions for reading variuos parts of gaussian cube files from char matrices
// checking geometry and preparing surfaces
void read_check_geom_surf_cube(SURFACE *surf, int s, ATOMS *xyz, int natoms, int centerofmass, double *centerofmasscoor, char **input)
{
int i, j, k, c, nactatom, tmpatomnum;
char *string, *token;
double tmpcoor[3];

if(input[0]==NULL || input[1]==NULL || input[2]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[2], "%d", &nactatom))!=1) {puts("\nProblems when reading xyz data...\n"); exit(1);}
if(nactatom<0) {nactatom = -nactatom; surf[s].sign = 2;}
 else surf[s].sign = 1;
if(nactatom<1) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if(nactatom!=natoms) {puts("\nGeometries from cube files differ...\nUse only cube files for one geometry of a molecule in one charmol input!\n"); exit(1);}
if(input[3]==NULL || input[4]==NULL || input[5]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
for(j=0,i=6;j<nactatom;j++,i++)
 {
 if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 if((c = sscanf(input[i], "%d%*f%lf%lf%lf", &tmpatomnum, &tmpcoor[0], &tmpcoor[1], &tmpcoor[2]))!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}
 if(centerofmass==1) {for(c=0;c<3;c++) tmpcoor[c] -= centerofmasscoor[c];}
 if(xyz[j].atomnum!=tmpatomnum) {puts("\nGeometries from cube files differ...\nUse only cube files for one geometry of a molecule in one charmol input!\n"); exit(1);}
 for(c=0;c<3;c++) {if(fabs(xyz[j].coor[c] - tmpcoor[c])>1E-5)
  {puts("\nGeometries from cube files differ...\nUse only cube files for one geometry of a molecule in one charmol input!\n"); exit(1);}}
 }
if(surf[s].sign==1)
 {
 if(surf[s].nsurf!=0) puts("Ignoring user-specified orbitals settings for cube file not containing orbitals...");
 surf[s].nsurf = 1;
 if((surf[s].surfs = (int *) realloc(surf[s].surfs, 1*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
 surf[s].surfs[0] = 0;
 surf[s].norb = 0;
 surf[s].orbs = NULL;
 }
 else
  {
  if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
  string = input[i];
  if((token = read_next_token(string))==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
  string = NULL;
  if((c = sscanf(token, "%d", &surf[s].norb))!=1) {puts("\nProblems when reading cube orbitals data...\n"); exit(1);}
  if(surf[s].norb<1) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
  if((surf[s].orbs = (int *) malloc(surf[s].norb*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
  for(j=0;j<surf[s].norb;j++)
   {
   // orbital specifications can be on more lines
   if((token = read_next_token(string))==NULL)
    {
    // go to next line
    i++;
    if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
    string = input[i];
    if((token = read_next_token(string))==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
    string = NULL;
    }
   if((c = sscanf(token, "%d", &surf[s].orbs[j]))!=1) {puts("\nProblems when reading cube orbitals data...\n"); exit(1);}
   }
  if(surf[s].nsurf>surf[s].norb) {puts("\nYou require depiction of more orbitals than found in cube file...\n"); exit(1);}
  if(surf[s].nsurf==0)
   {
   surf[s].nsurf = surf[s].norb;
   if((surf[s].surfs = (int *) realloc(surf[s].surfs, surf[s].nsurf*sizeof(int)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
   for(j=0;j<surf[s].nsurf;j++) surf[s].surfs[j] = j;
   }
   else
    {
    for(j=0;j<surf[s].nsurf;j++)
     {
     for(k=0;k<surf[s].norb;k++) if(surf[s].surfs[j]==surf[s].orbs[k]) {surf[s].surfs[j] = k; break;}
     if(k==surf[s].norb) {printf("\nOrbital %d not found in cube file...\n\n", surf[s].surfs[j]); exit(1);}
     }
    }
  }
}

// reading grid information
void read_grid_info_cube(int *gpoints, double *gmin, double *griddist, int centerofmass, double *centerofmasscoor, char **input)
{
int c;
double dval[3];

if(input[0]==NULL || input[1]==NULL || input[2]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[2], "%*d%lf%lf%lf", &gmin[0], &gmin[1], &gmin[2]))!=3) {puts("\nProblems when reading grid data...\n"); exit(1);}
if(centerofmass==1) {for(c=0;c<3;c++) gmin[c] -= centerofmasscoor[c];}
if(input[3]==NULL || input[4]==NULL || input[5]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[3], "%d%lf%lf%lf", &gpoints[0], &dval[0], &dval[1], &dval[2]))!=4) {puts("\nProblems when reading grid data...\n"); exit(1);}
if(fabs(dval[1])>9E-7 || fabs(dval[2])>9E-7) {puts("\nOnly orthogonal voxels are supported in charmol program...\n"); exit(1);}
griddist[0] = dval[0];
if((c = sscanf(input[4], "%d%lf%lf%lf", &gpoints[1], &dval[0], &dval[1], &dval[2]))!=4) {puts("\nProblems when reading grid data...\n"); exit(1);}
if(fabs(dval[0])>9E-7 || fabs(dval[2])>9E-7) {puts("\nOnly orthogonal voxels are supported in charmol program...\n"); exit(1);}
griddist[1] = dval[1];
if((c = sscanf(input[5], "%d%lf%lf%lf", &gpoints[2], &dval[0], &dval[1], &dval[2]))!=4) {puts("\nProblems when reading grid data...\n"); exit(1);}
if(fabs(dval[0])>9E-7 || fabs(dval[1])>9E-7) {puts("\nOnly orthogonal voxels are supported in charmol program...\n"); exit(1);}
griddist[2] = dval[2];
}

// checking geometry and grid settings
void read_check_geom_grid_cube(ATOMS *xyz, int natoms, int centerofmass, double *centerofmasscoor, int *gpoints, double *gmin, double *griddist, char **input)
{
int i, j, c, tmpnatoms, tmpatomnum, tmpgpoints;
double tmpcoor[3];

if(input[0]==NULL || input[1]==NULL || input[2]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[2], "%d%lf%lf%lf", &tmpnatoms, &tmpcoor[0], &tmpcoor[1], &tmpcoor[2]))!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}
if(centerofmass==1) {for(c=0;c<3;c++) tmpcoor[c] -= centerofmasscoor[c];}
if(fabs(gmin[0] - tmpcoor[0])>9E-7 || fabs(gmin[1] - tmpcoor[1])>9E-7 || fabs(gmin[2] - tmpcoor[2])>9E-7)
 {puts("\nGrid settings for potential and surface, which you want to make colormapped, must be exactly the same...\n"); exit(1);}
if(tmpnatoms<0) tmpnatoms = -tmpnatoms;
if(tmpnatoms<1) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if(tmpnatoms!=natoms) {puts("\nGeometries from cube files differ...\nUse only cube files for one geometry of a molecule in one charmol input!\n"); exit(1);}
if(input[3]==NULL || input[4]==NULL || input[5]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[3], "%d%lf%lf%lf", &tmpgpoints, &tmpcoor[0], &tmpcoor[1], &tmpcoor[2]))!=4) {puts("\nProblems when reading grid data...\n"); exit(1);}
if(fabs(tmpcoor[1])>9E-7 || fabs(tmpcoor[2])>9E-7) {puts("\nOnly orthogonal voxels are supported in charmol program...\n"); exit(1);}
if(tmpgpoints!=gpoints[0] || fabs(griddist[0] - tmpcoor[0])>9E-7)
 {puts("\nGrid settings for potential and surface, which you want to make colormapped, must be exactly the same...\n"); exit(1);}
if((c = sscanf(input[4], "%d%lf%lf%lf", &tmpgpoints, &tmpcoor[0], &tmpcoor[1], &tmpcoor[2]))!=4) {puts("\nProblems when reading grid data...\n"); exit(1);}
if(fabs(tmpcoor[0])>9E-7 || fabs(tmpcoor[2])>9E-7) {puts("\nOnly orthogonal voxels are supported in charmol program...\n"); exit(1);}
if(tmpgpoints!=gpoints[1] || fabs(griddist[1] - tmpcoor[1])>9E-7)
 {puts("\nGrid settings for potential and surface, which you want to make colormapped, must be exactly the same...\n"); exit(1);}
if((c = sscanf(input[5], "%d%lf%lf%lf", &tmpgpoints, &tmpcoor[0], &tmpcoor[1], &tmpcoor[2]))!=4) {puts("\nProblems when reading grid data...\n"); exit(1);}
if(fabs(tmpcoor[0])>9E-7 || fabs(tmpcoor[1])>9E-7) {puts("\nOnly orthogonal voxels are supported in charmol program...\n"); exit(1);}
if(tmpgpoints!=gpoints[2] || fabs(griddist[2] - tmpcoor[2])>9E-7)
 {puts("\nGrid settings for potential and surface, which you want to make colormapped, must be exactly the same...\n"); exit(1);}
for(j=0,i=6;j<tmpnatoms;j++,i++)
 {
 if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 if((c = sscanf(input[i], "%d%*f%lf%lf%lf", &tmpatomnum, &tmpcoor[0], &tmpcoor[1], &tmpcoor[2]))!=4) {puts("\nProblems when reading xyz data...\n"); exit(1);}
 if(centerofmass==1) {for(c=0;c<3;c++) tmpcoor[c] -= centerofmasscoor[c];}
 if(xyz[j].atomnum!=tmpatomnum) {puts("\nGeometries from cube files differ...\nUse only cube files for one geometry of a molecule in one charmol input!\n"); exit(1);}
 for(c=0;c<3;c++) {if(fabs(xyz[j].coor[c] - tmpcoor[c])>1E-5)
  {puts("\nGeometries from cube files differ...\nUse only cube files for one geometry of a molecule in one charmol input!\n"); exit(1);}}
 }
}

// reading one set of volumetric data - values in points of grid
void read_grid_values_oneset_cube(double *gridval, int totpoints, char **input)
{
int i, j, c, nat, orbs, norb, junk;
char *string, *token;

if(input[0]==NULL || input[1]==NULL || input[2]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[2], "%d", &nat))!=1) {puts("\nProblems when reading grid data...\n"); exit(1);}
if(nat<0) {nat = -nat; orbs = 1;}
 else orbs = 0;
if(input[3]==NULL || input[4]==NULL || input[5]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
for(j=0,i=6;j<nat;j++,i++) if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if(orbs==1)
 {
 // orbital specifications can be on more lines (we have to skip them all)
 if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%d", &norb))!=1) {puts("\nProblems when reading cube orbitals data...\n"); exit(1);}
 if(norb<1) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 for(j=0;j<norb;j++)
  {
  if((token = read_next_token(string))==NULL)
   {
   // go to next line
   i++;
   if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
   string = input[i];
   if((token = read_next_token(string))==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
   string = NULL;
   }
  if((c = sscanf(token, "%d", &junk))!=1) {puts("\nProblems when reading cube orbitals data...\n"); exit(1);}
  }
 }
 else i--; // go back to last line before volumetric data
j = 0;
while(j<totpoints)
 {
 i++; if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL)
  {
  // there can be one empty line between blocks of volumetric data in cube file
  i++; if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
  string = input[i];
  if((token = read_next_token(string))==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
  }
 string = NULL;
 if((c = sscanf(token, "%lf", &gridval[j]))!=1) {puts("\nProblems when reading grid data...\n"); exit(1);}
 j++;
 while((token = read_next_token(string))!=NULL && j<totpoints)
  {
  if((c = sscanf(token, "%lf", &gridval[j]))!=1) {puts("\nProblems when reading grid data...\n"); exit(1);}
  j++;
  }
 }
}

// reading more sets of volumetric data - values in points of grid
void read_grid_values_cube(double **gridval, int totpoints, int totsurfs, char **input)
{
int i, j, k, l, c, nat, orbs, norb, junk;
char *string, *token;

if(input[0]==NULL || input[1]==NULL || input[2]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[2], "%d", &nat))!=1) {puts("\nProblems when reading grid data...\n"); exit(1);}
if(nat<0) {nat = -nat; orbs = 1;}
 else orbs = 0;
if(input[3]==NULL || input[4]==NULL || input[5]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
for(j=0,i=6;j<nat;j++,i++) if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
if(orbs==1)
 {
 // orbital specifications can be on more lines (we have to skip them all)
 if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 string = NULL;
 if((c = sscanf(token, "%d", &norb))!=1) {puts("\nProblems when reading cube orbitals data...\n"); exit(1);}
 if(norb<1) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 for(j=0;j<norb;j++)
  {
  if((token = read_next_token(string))==NULL)
   {
   // go to next line
   i++;
   if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
   string = input[i];
   if((token = read_next_token(string))==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
   string = NULL;
   }
  if((c = sscanf(token, "%d", &junk))!=1) {puts("\nProblems when reading cube orbitals data...\n"); exit(1);}
  }
 }
 else i--; // go back to last line before volumetric data
j = 0; k = 0; l = 0;
while(j<(totpoints*totsurfs))
 {
 i++; if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
 string = input[i];
 if((token = read_next_token(string))==NULL)
  {
  // there can be one empty line between blocks of volumetric data in cube file
  i++; if(input[i]==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
  string = input[i];
  if((token = read_next_token(string))==NULL) {puts("\nCube file seems to be corrupted...\n"); exit(1);}
  }
 string = NULL;
 if((c = sscanf(token, "%lf", &gridval[k][l]))!=1) {puts("\nProblems when reading grid data...\n"); exit(1);}
 j++; k++; if(k==totsurfs) {k = 0; l++;}
 while((token = read_next_token(string))!=NULL && j<(totpoints*totsurfs))
  {
  if((c = sscanf(token, "%lf", &gridval[k][l]))!=1) {puts("\nProblems when reading grid data...\n"); exit(1);}
  j++; k++; if(k==totsurfs) {k = 0; l++;}
  }
 }
}

// checking signs of actual grid values for non-orbital volumetric data (it can be difference density)
void read_check_grid_value_signs_cube(SURFACE *surf, int s, double **gridval, int totpoints, int totsurfs)
{
int i, j, actsign;

if(surf[s].sign==1) // check done only for non-orbital data
 {
 actsign = surf[s].sign;
 for(i=0;i<totsurfs;i++)
  for(j=0;j<totpoints;j++)
   if(gridval[i][j] < 0.0) actsign = 2; // negative value found
 if(actsign==2) surf[s].sign = actsign;
 }
}

/*********************************************************************************************************************************************************************************/
// functions for reading vibrations from char matrices
// reading vibration from molden (displacements and FR-COORD xyz in atomic units)
void read_vibration_molden(TYCKA *vibtyc, KONUS *vibkon, ATOMS *xyz, int natoms, int centerofmass, double *centerofmasscoor, int vibnum, char **input)
{
int i, j, c, num;
static int firsttime = 1;
double tmpcoor[3];

for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "[FR-NORM-COORD]"))!=NULL || (strstr(input[i], "[Fr-norm-coord]"))!=NULL || (strstr(input[i], "[fr-norm-coord]"))!=NULL) break; i++;
if(input[i-1]==NULL) {puts("\nNo vibrations found in molden file...\n"); exit(1);}
if((strstr(input[i], "vibration"))==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
if((c = sscanf(input[i], "%*s%d", &num))!=1) {puts("\nProblems when reading vibration data...\n"); exit(1);} i++;
while(num!=vibnum)
 {
 for(j=0;j<natoms;j++) {if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);} else i++;}
 if(input[i]==NULL) {printf("\nVibration %d not found in molden file...\n\n", vibnum); exit(1);}
 if((strstr(input[i], "vibration"))==NULL) {printf("\nVibration %d not found in molden file...\n\n", vibnum); exit(1);}
 if((c = sscanf(input[i], "%*s%d", &num))!=1) {puts("\nProblems when reading vibration data...\n"); exit(1);} i++;
 }
for(j=0;j<natoms;j++)
 {
 if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
 if((c = sscanf(input[i], "%lf%lf%lf", &vibkon[j].end[0], &vibkon[j].end[1], &vibkon[j].end[2]))!=3) {puts("\nProblems when reading vibration data...\n"); exit(1);}
 i++;
 }
if(firsttime==1)
 {
 for(j=0;j<natoms;j++)
  {
  for(c=0;c<3;c++) vibtyc[j].start[c] = xyz[j].coor[c];
  vibtyc[j].ordid = xyz[j].ordid;
  }
 for(i=0;input[i]!=NULL;i++) if((strstr(input[i], "[FR-COORD]"))!=NULL || (strstr(input[i], "[Fr-coord]"))!=NULL || (strstr(input[i], "[fr-coord]"))!=NULL) break; i++;
 if(input[i-1]!=NULL)
  {
  num = 0;
  for(j=0;j<natoms;j++)
   {
   if(input[i]==NULL) {puts("\nMolden file seems to be corrupted...\n"); exit(1);}
   if((c = sscanf(input[i], "%*s%lf%lf%lf", &tmpcoor[0], &tmpcoor[1], &tmpcoor[2]))!=3) {puts("\nProblems when reading vibration data...\n"); exit(1);}
   if(centerofmass==1) {for(c=0;c<3;c++) tmpcoor[c] -= centerofmasscoor[c];}
   for(c=0;c<3;c++) {if(fabs(vibtyc[j].start[c] - tmpcoor[c])>1E-3) num = 1;}
   i++;
   }
  if(num==1) puts("ATOMS and FR-COORD geometries from molden differ. ATOMS geometry is used for depiction of vibrations, so be careful about the result...");
  }
 firsttime = 0;
 }
}

