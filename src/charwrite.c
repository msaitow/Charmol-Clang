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
#include "charwrite.h"

/*********************************************************************************************************************************************************************************/
// functions for running system commands (external programs)
// function for running basic system command repaired for special (shell-interpreted) characters
void write_run_system_command(char *inputcommand)
{
int i, pos;
char command[600];

strcpy(command, inputcommand);

i = 0;
while(command[i]!='\0')
 {
 while(command[i]!='\'' && command[i]!='\"' && command[i]!='\\' && command[i]!='\0') i++;
 if(command[i]=='\'' || command[i]=='\"' || command[i]=='\\')
  {
  pos = i;
  while(command[i]!='\0') i++;
  while(i >= pos) {command[i+1] = command[i]; i--;}
  command[pos] = '\\';
  i = pos + 2;
  }
 }

system(command);
}

/*********************************************************************************************************************************************************************************/
// functions helping with determining stuff necessary for writing output
// function for determining vrml rotation (reference vector lies on y axis: 0 1 0)
void write_determine_vrml_rotation(double *start, double *end, double *rot)
{
int i;
double vec[3], veclen, scalar;

veclen = 0.0;
for(i=0;i<3;i++)
 {
 vec[i] = end[i] - start[i];
 veclen += vec[i]*vec[i];
 }
veclen = sqrt(veclen);
if(veclen>1E-6) scalar = vec[1]/veclen;
 else scalar = 0.0;
if(fabs(scalar)>(1.0 - 1E-6))
 {
 rot[0] = 1.0; rot[1] = 0.0; rot[2] = 0.0;
 if((end[1] - start[1])>0.0) rot[3] = 0.0;
  else rot[3] = M_PI;
 }
 else
  {
  veclen = sqrt(vec[0]*vec[0] + vec[2]*vec[2]);
  rot[0] = vec[2]/veclen; rot[1] = 0.0; rot[2] = -vec[0]/veclen; rot[3] = acos(scalar);
  }
}

/*********************************************************************************************************************************************************************************/
// functions for writing various types of 3D objects
// used outformat is: 1 - povray, 2 - vrml
// function for writing sphere
FILE *write_sphere(double *center, double rad, double *rgb, int finish, FILE *outfile, int outformat)
{

switch(outformat)
 {
 case 1:
  fprintf(outfile, "sphere\n {\n");
  fprintf(outfile, " <%15f,%15f,%15f>,%15f\n", center[0], center[1], center[2], rad);
  fprintf(outfile, " texture { pigment { rgb <%7.3f%7.3f%7.3f> }", rgb[0], rgb[1], rgb[2]);
  if(finish==1) fprintf(outfile, " finish { ambient BSAMBI diffuse BSDIFF specular BSSPEC } }\n }\n");
   else fprintf(outfile, " }\n }\n");
  break;
 case 2:
  fprintf(outfile, "Transform\n {\n");
  fprintf(outfile, " translation %15f%15f%15f\n", center[0], center[1], center[2]);
  fprintf(outfile, " children\n  [\n  Shape\n   {\n");
  fprintf(outfile, "   appearance Appearance { material charmolMaterial { diffuseColor %.3f %.3f %.3f } }\n", rgb[0], rgb[1], rgb[2]);
  fprintf(outfile, "   geometry Sphere { radius %f }\n   }\n  ]\n }\n", rad);
  break;
 }

return outfile;
}

// function for writing cylinder
FILE *write_cylinder(double *start, double *end, double *rot, double rad, double *rgb, int finish, FILE *outfile, int outformat)
{
int i;
double center[3], length;

switch(outformat)
 {
 case 1:
  fprintf(outfile, "cylinder\n {\n");
  fprintf(outfile, " <%15f,%15f,%15f>, <%15f,%15f,%15f>,%15f\n", start[0], start[1], start[2], end[0], end[1], end[2], rad);
  fprintf(outfile, " texture { pigment { rgb <%7.3f%7.3f%7.3f> }", rgb[0], rgb[1], rgb[2]);
  if(finish==1) fprintf(outfile, " finish { ambient BSAMBI diffuse BSDIFF specular BSSPEC } }\n }\n");
   else fprintf(outfile, " }\n }\n");
  break;
 case 2:
  length = 0.0;
  for(i=0;i<3;i++)
   {
   center[i] = start[i] + (end[i] - start[i])/2.0;
   length += (end[i] - start[i])*(end[i] - start[i]);
   }
  length = sqrt(length);
  fprintf(outfile, "Transform\n {\n");
  fprintf(outfile, " translation %15f%15f%15f\n", center[0], center[1], center[2]);
  fprintf(outfile, " rotation %15f%15f%15f%15f\n", rot[0], rot[1], rot[2], rot[3]);
  fprintf(outfile, " children\n  [\n  Shape\n   {\n");
  fprintf(outfile, "   appearance Appearance { material charmolMaterial { diffuseColor %.3f %.3f %.3f } }\n", rgb[0], rgb[1], rgb[2]);
  fprintf(outfile, "   geometry Cylinder { radius %f height %f top FALSE bottom FALSE }\n   }\n  ]\n }\n", rad, length);
  break;
 }

return outfile;
}

// function for writing cylinder with top and bottom
FILE *write_cylinder_topbot(double *start, double *end, double *rot, double rad, double *rgb, int finish, FILE *outfile, int outformat)
{
int i;
double center[3], length;

switch(outformat)
 {
 case 1:
  fprintf(outfile, "cylinder\n {\n");
  fprintf(outfile, " <%15f,%15f,%15f>, <%15f,%15f,%15f>,%15f\n", start[0], start[1], start[2], end[0], end[1], end[2], rad);
  fprintf(outfile, " texture { pigment { rgb <%7.3f%7.3f%7.3f> }", rgb[0], rgb[1], rgb[2]);
  if(finish==1) fprintf(outfile, " finish { ambient BSAMBI diffuse BSDIFF specular BSSPEC } }\n }\n");
   else fprintf(outfile, " }\n }\n");
  break;
 case 2:
  length = 0.0;
  for(i=0;i<3;i++)
   {
   center[i] = start[i] + (end[i] - start[i])/2.0;
   length += (end[i] - start[i])*(end[i] - start[i]);
   }
  length = sqrt(length);
  fprintf(outfile, "Transform\n {\n");
  fprintf(outfile, " translation %15f%15f%15f\n", center[0], center[1], center[2]);
  fprintf(outfile, " rotation %15f%15f%15f%15f\n", rot[0], rot[1], rot[2], rot[3]);
  fprintf(outfile, " children\n  [\n  Shape\n   {\n");
  fprintf(outfile, "   appearance Appearance { material charmolMaterial { diffuseColor %.3f %.3f %.3f } }\n", rgb[0], rgb[1], rgb[2]);
  fprintf(outfile, "   geometry Cylinder { radius %f height %f }\n   }\n  ]\n }\n", rad, length);
  break;
 }

return outfile;
}

// function for writing cone
FILE *write_cone(double *start, double *end, double *rot, double rad, double *rgb, int finish, FILE *outfile, int outformat)
{
int i;
double center[3];

switch(outformat)
 {
 case 1:
  fprintf(outfile, "cone\n {\n");
  fprintf(outfile, " <%15f,%15f,%15f>,%15f, <%15f,%15f,%15f>,%15f\n", start[0], start[1], start[2], rad, end[0], end[1], end[2], 0.0);
  fprintf(outfile, " texture { pigment { rgb <%7.3f%7.3f%7.3f> }", rgb[0], rgb[1], rgb[2]);
  if(finish==1) fprintf(outfile, " finish { ambient BSAMBI diffuse BSDIFF specular BSSPEC } }\n }\n");
   else fprintf(outfile, " }\n }\n");
  break;
 case 2:
  for(i=0;i<3;i++) center[i] = start[i] + (end[i] - start[i])/2.0;
  fprintf(outfile, "Transform\n {\n");
  fprintf(outfile, " translation %15f%15f%15f\n", center[0], center[1], center[2]);
  fprintf(outfile, " rotation %15f%15f%15f%15f\n", rot[0], rot[1], rot[2], rot[3]);
  fprintf(outfile, " children\n  [\n  Shape\n   {\n");
  fprintf(outfile, "   appearance Appearance { material charmolMaterial { diffuseColor %.3f %.3f %.3f } }\n", rgb[0], rgb[1], rgb[2]);
  fprintf(outfile, "   geometry Cone { bottomRadius %f height %f }\n   }\n  ]\n }\n", rad, rad);
  break;
 }

return outfile;
}

/*********************************************************************************************************************************************************************************/
// main functions for writing various parts of the output
// used outformat is: 1 - povray, 2 - vrml
// writing basic settings
FILE *write_basic_settings(ATOMS *xyz, int natoms, TYCKA *tyc, KONUS *kon, int narrows, int camspec, int useangle, double *camera,
                           double convertxyzval, int scenerot, double scenerotangle, int camrot, double *camerarot, int camrq,
                           double *camerarq, int camsky, double *camerasky, int camzoom, double camerazoom, int handedness,
                           int finish, double *finset, FILE *outfile, int outformat)
{
int i, j;
double maxcoor, tmpcamera[3], tmpcamerasky[3], tmpcamerarq[4], rq[4][4], rqlen, tmpval[3];

switch(outformat)
 {
 case 1:
  if(camspec==0)
   {
   maxcoor = 0.0;
   for(i=0;i<natoms;i++)
    {
    if(fabs(xyz[i].coor[0])>maxcoor) maxcoor = fabs(xyz[i].coor[0]);
    if(fabs(xyz[i].coor[1])>maxcoor) maxcoor = fabs(xyz[i].coor[1]);
    }
   for(i=0;i<narrows;i++)
    {
    if(fabs(tyc[i].start[0])>maxcoor) maxcoor = fabs(tyc[i].start[0]);
    if(fabs(tyc[i].start[1])>maxcoor) maxcoor = fabs(tyc[i].start[1]);
    if(fabs(kon[i].end[0])>maxcoor) maxcoor = fabs(kon[i].end[0]);
    if(fabs(kon[i].end[1])>maxcoor) maxcoor = fabs(kon[i].end[1]);
    }
   if(maxcoor<0.001) maxcoor = 1.0;
   if(useangle==1) maxcoor *= 2.2;
    else maxcoor *= 0.45;
   tmpcamera[0] = 0.0; tmpcamera[1] = 0.0;
   tmpcamera[2] = maxcoor/tan(M_PI*7.5/180.0);
   }
   else
    {
    for(i=0;i<3;i++) tmpcamera[i] = camera[i]*convertxyzval;
    }
  if(camzoom==1) {for(i=0;i<3;i++) tmpcamera[i] *= camerazoom;}
  if(camsky==1) {for(i=0;i<3;i++) tmpcamerasky[i] = camerasky[i];}
   else {tmpcamerasky[0] = 0.0; tmpcamerasky[1] = 1.0; tmpcamerasky[2] = 0.0;}
  if(camrq==1)
   {
   rqlen = 0.0; for(i=0;i<4;i++) rqlen += camerarq[i]*camerarq[i]; rqlen = sqrt(rqlen);
   for(i=0;i<4;i++) tmpcamerarq[i] = camerarq[i]/rqlen;
   for(i=0;i<4;i++) {for(j=0;j<4;j++) rq[i][j] = tmpcamerarq[i]*tmpcamerarq[j];}
   tmpval[0] = (rq[0][0] + rq[1][1] - rq[2][2] - rq[3][3]) * tmpcamera[0] +
               (2.0*rq[1][2] + 2.0*rq[0][3])               * tmpcamera[1] +
               (2.0*rq[1][3] - 2.0*rq[0][2])               * tmpcamera[2];
   tmpval[1] = (2.0*rq[1][2] - 2.0*rq[0][3])               * tmpcamera[0] +
               (rq[0][0] - rq[1][1] + rq[2][2] - rq[3][3]) * tmpcamera[1] +
               (2.0*rq[2][3] + 2.0*rq[0][1])               * tmpcamera[2];
   tmpval[2] = (2.0*rq[1][3] + 2.0*rq[0][2])               * tmpcamera[0] +
               (2.0*rq[2][3] - 2.0*rq[0][1])               * tmpcamera[1] +
               (rq[0][0] - rq[1][1] - rq[2][2] + rq[3][3]) * tmpcamera[2];
   for(i=0;i<3;i++) tmpcamera[i] = tmpval[i];
   tmpval[0] = (rq[0][0] + rq[1][1] - rq[2][2] - rq[3][3]) * tmpcamerasky[0] +
               (2.0*rq[1][2] + 2.0*rq[0][3])               * tmpcamerasky[1] +
               (2.0*rq[1][3] - 2.0*rq[0][2])               * tmpcamerasky[2];
   tmpval[1] = (2.0*rq[1][2] - 2.0*rq[0][3])               * tmpcamerasky[0] +
               (rq[0][0] - rq[1][1] + rq[2][2] - rq[3][3]) * tmpcamerasky[1] +
               (2.0*rq[2][3] + 2.0*rq[0][1])               * tmpcamerasky[2];
   tmpval[2] = (2.0*rq[1][3] + 2.0*rq[0][2])               * tmpcamerasky[0] +
               (2.0*rq[2][3] - 2.0*rq[0][1])               * tmpcamerasky[1] +
               (rq[0][0] - rq[1][1] - rq[2][2] + rq[3][3]) * tmpcamerasky[2];
   for(i=0;i<3;i++) tmpcamerasky[i] = tmpval[i];
   }
  fprintf(outfile, "// povray generated by CHARMOL program\n// convertxyzval %f \n\n// basic settings\n", convertxyzval);
  if(scenerot==1) fprintf(outfile, "#include \"transforms.inc\"\n");
  fprintf(outfile, "background { color rgb <1.0, 1.0, 1.0> }\n");
  fprintf(outfile, "camera\n {\n location <%.3f, %.3f, %.3f>\n look_at <0.0, 0.0, 0.0>\n", tmpcamera[0], tmpcamera[1], tmpcamera[2]);
  if(camrot==1) fprintf(outfile, " rotate <%.3f, %.3f, %.3f>\n", camerarot[0], camerarot[1], camerarot[2]);
  if(useangle==1) fprintf(outfile, " angle 15\n");
  if(handedness==-1) fprintf(outfile, " right <%.2f, 0, 0>\n", 1.33*handedness);
  if(camsky==1 || camrq==1) fprintf(outfile, " sky <%.3f, %.3f, %.3f>\n", tmpcamerasky[0], tmpcamerasky[1], tmpcamerasky[2]);
  if(scenerot==1) fprintf(outfile, " Axis_Rotate_Trans(<%.3f, %.3f, %.3f>, %.3f)\n", tmpcamera[0], tmpcamera[1], tmpcamera[2], (-1)*handedness*scenerotangle);
  fprintf(outfile, " }\n");
  fprintf(outfile, "light_source { <%.3f, %.3f, %.3f> color rgb <1, 1, 1>", tmpcamera[0], tmpcamera[1], tmpcamera[2]);
  if(camrot==1) fprintf(outfile, " rotate <%.3f, %.3f, %.3f> }\n", camerarot[0], camerarot[1], camerarot[2]);
   else fprintf(outfile, " }\n");
  if(finish==1) fprintf(outfile, "#declare BSAMBI = %.3f;\n#declare BSDIFF = %.3f;\n#declare BSSPEC = %.3f;\n\n", finset[0], finset[1], finset[2]);
   else fprintf(outfile, "\n");
  break;
 case 2:
  fprintf(outfile, "#VRML V2.0 utf8\n");
  fprintf(outfile, "# vrml generated by CHARMOL program\n\n# basic settings\n");
  fprintf(outfile, "NavigationInfo { type \"EXAMINE\" }\n");
  fprintf(outfile, "Background { skyColor [ 1 1 1 ] }\n");
  fprintf(outfile, "PROTO charmolMaterial\n [\n");
  fprintf(outfile, " exposedField SFFloat ambientIntensity %.3f\n", finset[3]);
  fprintf(outfile, " exposedField SFColor diffuseColor     0.8 0.8 0.8\n");
  fprintf(outfile, " exposedField SFColor emissiveColor    0.0 0.0 0.0\n");
  fprintf(outfile, " exposedField SFFloat shininess        %.3f\n", finset[4]);
  if(finish==1) fprintf(outfile, " exposedField SFColor specularColor    1.0 1.0 1.0\n");
   else fprintf(outfile, " exposedField SFColor specularColor    0.0 0.0 0.0\n");
  fprintf(outfile, " exposedField SFFloat transparency     0.0\n");
  fprintf(outfile, " ]\n {\n Material\n  {\n");
  fprintf(outfile, "  ambientIntensity IS ambientIntensity\n  diffuseColor     IS diffuseColor\n");
  fprintf(outfile, "  specularColor    IS specularColor\n  emissiveColor    IS emissiveColor\n");
  fprintf(outfile, "  shininess        IS shininess\n  transparency     IS transparency\n  }\n }\n\n");
  break;
 }

return outfile;
}

// writing xyz
FILE *write_xyz(ATOMS *xyz, int natoms, int nmolecules, int hideatoms, int atomcolor, int bondcolor, int finish, double *atomradiusscale,
                double *lightatomradiusscale, double *bondradius, double **defatomcolor, double **defbondcolor, FILE *outfile, int outformat)
{
int i, j, k;
double rad1, rad2, dist, end[3], rot[4], rgb[3], bondrad;

for(i=0;i<natoms;i++)
 {
 // making atomic sphere
 switch(outformat)
  {
  case 1: fprintf(outfile, "// atom %d\n", i+1); break;
  case 2: fprintf(outfile, "# atom %d\n", i+1); break;
  }
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
 outfile = write_sphere(xyz[i].coor, rad1, rgb, finish, outfile, outformat);
 // making its halfbonds
 for(j=0;xyz[i].cone[j]!=0;j++)
  {
  write_determine_vrml_rotation(xyz[i].coor, xyz[xyz[i].cone[j]-1].coor, rot);
  switch(outformat)
   {
   case 1: fprintf(outfile, "// halfbond %d - %d\n", i+1, xyz[i].cone[j]); break;
   case 2: fprintf(outfile, "# halfbond %d - %d\n", i+1, xyz[i].cone[j]); break;
   }
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
  outfile = write_cylinder(xyz[i].coor, end, rot, bondrad, rgb, finish, outfile, outformat);
  }
 fprintf(outfile, "\n");
 }

return outfile;
}

// writing gauges
FILE *write_gauges(ATOMS *xyz, int natoms, GAUGE *gauge, int ngauges, int hideatoms, int finish, double *atomradiusscale,
                   double *lightatomradiusscale, double *bondradius, FILE *outfile, int outformat)
{
int i, j, l, c, a1, a2, extradash;
double rad1, rad2, dist, dash, ball, space, rest, start[3], end[3], rot[4];

for(i=0;i<ngauges;i++)
 {
 a1 = gauge[i].a1 - 1;
 a2 = gauge[i].a2 - 1;
 write_determine_vrml_rotation(xyz[a1].coor, xyz[a2].coor, rot);
 switch(outformat)
  {
  case 1: fprintf(outfile, "// gauge %d - %d\n", a1+1, a2+1); break;
  case 2: fprintf(outfile, "# gauge %d - %d\n", a1+1, a2+1); break;
  }
 if(gauge[i].style!=4)
  {
  if(hideatoms==1) {rad1 = bondradius[xyz[a1].ordid]; rad2 = bondradius[xyz[a2].ordid];}
   else
    {
    if(xyz[a1].rad<(0.5*angstrom_to_bohr))
     {if(xyz[a1].rad*xyz[a1].scale*atomradiusscale[xyz[a1].ordid]*lightatomradiusscale[xyz[a1].ordid] < bondradius[xyz[a1].ordid]) rad1 = bondradius[xyz[a1].ordid];
       else rad1 = xyz[a1].rad*xyz[a1].scale*atomradiusscale[xyz[a1].ordid]*lightatomradiusscale[xyz[a1].ordid];}
     else
      {if(xyz[a1].rad*xyz[a1].scale*atomradiusscale[xyz[a1].ordid] < bondradius[xyz[a1].ordid]) rad1 = bondradius[xyz[a1].ordid];
       else rad1 = xyz[a1].rad*xyz[a1].scale*atomradiusscale[xyz[a1].ordid];}
    if(xyz[a2].rad<(0.5*angstrom_to_bohr))
     {if(xyz[a2].rad*xyz[a2].scale*atomradiusscale[xyz[a2].ordid]*lightatomradiusscale[xyz[a2].ordid] < bondradius[xyz[a2].ordid]) rad2 = bondradius[xyz[a2].ordid];
       else rad2 = xyz[a2].rad*xyz[a2].scale*atomradiusscale[xyz[a2].ordid]*lightatomradiusscale[xyz[a2].ordid];}
     else
      {if(xyz[a2].rad*xyz[a2].scale*atomradiusscale[xyz[a2].ordid] < bondradius[xyz[a2].ordid]) rad2 = bondradius[xyz[a2].ordid];
       else rad2 = xyz[a2].rad*xyz[a2].scale*atomradiusscale[xyz[a2].ordid];}
    }
  dist = 0.0; for(j=0;j<3;j++) dist += (xyz[a2].coor[j] - xyz[a1].coor[j])*(xyz[a2].coor[j] - xyz[a1].coor[j]); dist = sqrt(dist);
  }
 switch(gauge[i].style)
  {
  case 1: // dashed
   dash = 6.0*gauge[i].rad; space = 2.0*gauge[i].rad;
   c = (int) ((dist - rad1 - rad2)/(dash + space));
   rest = dist - rad1 - rad2 - (double) c*(dash + space);
   space += (double) rest/c;
   if(c==0) {puts("\nUnpossible to make dashed gauge. Check its thickness, which you require (make it less thick)...\n"); exit(1);}
   for(j=0;j<3;j++) {start[j] = xyz[a1].coor[j] + ((rad1 + space/2.0)/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]); end[j] = start[j] + (dash/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]);}
   for(l=0;l<c;l++)
    {
    outfile = write_cylinder_topbot(start, end, rot, gauge[i].rad, gauge[i].rgb, finish, outfile, outformat);
    for(j=0;j<3;j++) {start[j] = end[j] + (space/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]); end[j] = start[j] + (dash/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]);}
    }
   break;
  case 2: // dotted
   ball = 2.0*gauge[i].rad; space = 2.0*gauge[i].rad;
   c = (int) ((dist - rad1 - rad2)/(ball + space));
   rest = dist - rad1 - rad2 - (double) c*(ball + space);
   space += (double) rest/c;
   if(c==0) {puts("\nUnpossible to make dotted gauge. Check its thickness, which you require (make it less thick)...\n"); exit(1);}
   for(j=0;j<3;j++) start[j] = xyz[a1].coor[j] + ((rad1 + space/2.0 + ball/2.0)/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]);
   for(l=0;l<c;l++)
    {
    outfile = write_sphere(start, gauge[i].rad, gauge[i].rgb, finish, outfile, outformat);
    for(j=0;j<3;j++) start[j] += ((space + ball)/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]);
    }
   break;
  case 3: // dotdashed
   dash = 6.0*gauge[i].rad; ball = 2.0*gauge[i].rad; space = 1.5*gauge[i].rad;
   c = (int) ((dist - rad1 - rad2)/(dash + ball + 2.0*space));
   rest = dist - rad1 - rad2 - (double) c*(dash + ball + 2.0*space);
   if(rest < (dash + space)) {space += (double) rest/(2.0*c); extradash = 0;}
    else {rest -= dash + space; space += (double) rest/(2.0*c + 1.0); extradash = 1;}
   if(c==0) {puts("\nUnpossible to make dashed gauge. Check its thickness, which you require (make it less thick)...\n"); exit(1);}
   for(j=0;j<3;j++) {start[j] = xyz[a1].coor[j] + ((rad1 + space/2.0)/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]); end[j] = start[j] + (dash/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]);}
   for(l=0;l<c;l++)
    {
    outfile = write_cylinder_topbot(start, end, rot, gauge[i].rad, gauge[i].rgb, finish, outfile, outformat);
    for(j=0;j<3;j++) start[j] = end[j] + ((space + ball/2.0)/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]);
    outfile = write_sphere(start, gauge[i].rad, gauge[i].rgb, finish, outfile, outformat);
    for(j=0;j<3;j++) {start[j] += ((ball/2.0 + space)/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]); end[j] = start[j] + (dash/dist)*(xyz[a2].coor[j] - xyz[a1].coor[j]);}
    }
   if(extradash==1)
    {
    outfile = write_cylinder_topbot(start, end, rot, gauge[i].rad, gauge[i].rgb, finish, outfile, outformat);
    }
   break;
  case 4: // solid
   outfile = write_cylinder(xyz[a1].coor, xyz[a2].coor, rot, gauge[i].rad, gauge[i].rgb, finish, outfile, outformat);
   break;
  }
 fprintf(outfile, "\n");
 }

return outfile;
}

// writing single arrow
FILE *write_single_arrow(TYCKA tyc, KONUS kon, int finish, FILE *outfile, int outformat)
{
double rot[4];

write_determine_vrml_rotation(tyc.start, kon.end, rot);
switch(outformat)
 {
 case 1: fprintf(outfile, "// arrow\n"); break;
 case 2: fprintf(outfile, "# arrow\n"); break;
 }
outfile = write_cylinder(tyc.start, tyc.end, rot, tyc.rad, tyc.rgb, finish, outfile, outformat);
outfile = write_cone(tyc.end, kon.end, rot, kon.rad, tyc.rgb, finish, outfile, outformat);
fprintf(outfile, "\n");

return outfile;
}

// writing more arrows
FILE *write_arrows(TYCKA *tyc, KONUS *kon, int narrows, int finish, FILE *outfile, int outformat)
{
int i;

for(i=0;i<narrows;i++) outfile = write_single_arrow(tyc[i], kon[i], finish, outfile, outformat);

return outfile;
}

// writing solid mesh
FILE *write_mesh_solid(MESH mesh, double *rgb, double transparency, double isoval, int finish, FILE *outfile, int outformat)
{
int i, j;

if(mesh.ntriangles==0) return outfile;

switch(outformat)
 {
 case 1:
  fprintf(outfile, "// mesh\nmesh2\n {\n");
  fprintf(outfile, " vertex_vectors\n  {\n  %d,\n", 3*mesh.ntriangles);
  for(i=0;i<mesh.ntriangles;i++)
   {
   fprintf(outfile, "  <%16f,%16f,%16f>,\n  <%16f,%16f,%16f>,\n  <%16f,%16f,%16f>,\n",
           mesh.corners[i][0], mesh.corners[i][1], mesh.corners[i][2],
           mesh.corners[i][3], mesh.corners[i][4], mesh.corners[i][5],
           mesh.corners[i][6], mesh.corners[i][7], mesh.corners[i][8]);
   }
  fprintf(outfile, "  }\n");
  fprintf(outfile, " normal_vectors\n  {\n  %d,\n", 3*mesh.ntriangles);
  for(i=0;i<mesh.ntriangles;i++)
   {
   fprintf(outfile, "  <%16f,%16f,%16f>,\n  <%16f,%16f,%16f>,\n  <%16f,%16f,%16f>,\n",
           mesh.normals[i][0], mesh.normals[i][1], mesh.normals[i][2],
           mesh.normals[i][3], mesh.normals[i][4], mesh.normals[i][5],
           mesh.normals[i][6], mesh.normals[i][7], mesh.normals[i][8]);
   }
  fprintf(outfile, "  }\n");
  if(mesh.colors!=NULL)
   {
   fprintf(outfile, " texture_list\n  {\n  %d,\n", 3*mesh.ntriangles);
   for(i=0;i<mesh.ntriangles;i++)
    {
    for(j=0;j<3;j++)
     {
     fprintf(outfile, "  texture { pigment { rgb <%7.3f%7.3f%7.3f>", mesh.colors[i][3*j], mesh.colors[i][3*j+1], mesh.colors[i][3*j+2]);
     if(transparency>1E-3) fprintf(outfile, " transmit %.3f }", transparency);
      else fprintf(outfile, " }");
     if(finish==1) fprintf(outfile, " finish { ambient BSAMBI diffuse BSDIFF specular BSSPEC } }\n");
      else fprintf(outfile, " }\n");
     }
    }
   fprintf(outfile, "  }\n");
   fprintf(outfile, " face_indices\n  {\n  %d,\n", mesh.ntriangles);
   for(i=0,j=0;i<mesh.ntriangles;i++)
    {
    fprintf(outfile, "  <%d,%d,%d>,%d,%d,%d,\n", j, j+1, j+2, j, j+1, j+2);
    j += 3;
    }
   fprintf(outfile, "  }\n");
   }
   else
    {
    fprintf(outfile, " face_indices\n  {\n  %d,\n", mesh.ntriangles);
    for(i=0,j=0;i<mesh.ntriangles;i++)
     {
     fprintf(outfile, "  <%d,%d,%d>,\n", j, j+1, j+2);
     j += 3;
     }
    fprintf(outfile, "  }\n");
    fprintf(outfile, " texture { pigment { rgb <%7.3f%7.3f%7.3f>", rgb[0], rgb[1], rgb[2]);
    if(transparency>1E-3) fprintf(outfile, " transmit %.3f }", transparency);
     else fprintf(outfile, " }");
    if(finish==1) fprintf(outfile, " finish { ambient BSAMBI diffuse BSDIFF specular BSSPEC } }\n");
     else fprintf(outfile, " }\n");
    }
  fprintf(outfile, " }\n\n");
  break;
 case 2:
  fprintf(outfile, "# mesh\nShape\n {\n");
  fprintf(outfile, " appearance Appearance { material charmolMaterial {");
  if(mesh.colors==NULL) fprintf(outfile, " diffuseColor %.3f %.3f %.3f", rgb[0], rgb[1], rgb[2]);
  if(transparency>1E-3) fprintf(outfile, " transparency %.3f } }\n", transparency);
   else fprintf(outfile, " } }\n");
  fprintf(outfile, " geometry IndexedFaceSet\n  {\n");
  fprintf(outfile, "  coord Coordinate\n   {\n   point\n    [\n");
  for(i=0;i<mesh.ntriangles;i++)
   {
   fprintf(outfile, "%16f%16f%16f,\n%16f%16f%16f,\n%16f%16f%16f,\n",
           mesh.corners[i][0], mesh.corners[i][1], mesh.corners[i][2],
           mesh.corners[i][3], mesh.corners[i][4], mesh.corners[i][5],
           mesh.corners[i][6], mesh.corners[i][7], mesh.corners[i][8]);
   }
  fprintf(outfile, "    ]\n   }\n");
  fprintf(outfile, "  normal Normal\n   {\n   vector\n    [\n");
  for(i=0;i<mesh.ntriangles;i++)
   {
   if(isoval>0.0) fprintf(outfile, "%16f%16f%16f,\n%16f%16f%16f,\n%16f%16f%16f,\n",
                          -mesh.normals[i][0], -mesh.normals[i][1], -mesh.normals[i][2],
                          -mesh.normals[i][3], -mesh.normals[i][4], -mesh.normals[i][5],
                          -mesh.normals[i][6], -mesh.normals[i][7], -mesh.normals[i][8]);
    else fprintf(outfile, "%16f%16f%16f,\n%16f%16f%16f,\n%16f%16f%16f,\n",
                 mesh.normals[i][0], mesh.normals[i][1], mesh.normals[i][2],
                 mesh.normals[i][3], mesh.normals[i][4], mesh.normals[i][5],
                 mesh.normals[i][6], mesh.normals[i][7], mesh.normals[i][8]);
   }
  fprintf(outfile, "    ]\n   }\n");
  fprintf(outfile, "  coordIndex\n   [\n");
  for(i=0,j=0;i<mesh.ntriangles;i++)
   {
   if(isoval>0.0) fprintf(outfile, "%16d%16d%16d%16d,\n", j+2, j+1, j, -1);
    else fprintf(outfile, "%16d%16d%16d%16d,\n", j, j+1, j+2, -1);
   j += 3;
   }
  fprintf(outfile, "   ]\n");
  if(mesh.colors!=NULL)
   {
   fprintf(outfile, "  color Color\n   {\n   color\n    [\n");
   for(i=0;i<mesh.ntriangles;i++)
    {
    fprintf(outfile, "%16f%16f%16f,\n%16f%16f%16f,\n%16f%16f%16f,\n",
            mesh.colors[i][0], mesh.colors[i][1], mesh.colors[i][2],
            mesh.colors[i][3], mesh.colors[i][4], mesh.colors[i][5],
            mesh.colors[i][6], mesh.colors[i][7], mesh.colors[i][8]);
    }
   fprintf(outfile, "    ]\n   }\n");
   }
  fprintf(outfile, "  }\n }\n\n");
  break;
 }

return outfile;
}

// writing wireframe mesh
FILE *write_mesh_wireframe(MESH mesh, double *rgb, double transparency, double isoval, int finish, FILE *outfile, int outformat)
{
int i, j;
double *r1, *r2, *l1, *l2, rad;

if(mesh.ntriangles==0) return outfile;

switch(outformat)
 {
 case 1:
  fprintf(outfile, "// mesh\n");
  rad = 0.01;
  if(mesh.colors!=NULL)
   {
   for(i=0;i<mesh.ntriangles;i++)
    {
    for(j=0;j<3;j++)
     {
     switch(j)
      {
      case 0: r1 = &mesh.corners[i][0]; r2 = &mesh.corners[i][3]; l1 = &mesh.colors[i][0]; l2 = &mesh.colors[i][3]; break;
      case 1: r1 = &mesh.corners[i][0]; r2 = &mesh.corners[i][6]; l1 = &mesh.colors[i][0]; l2 = &mesh.colors[i][6]; break;
      case 2: r1 = &mesh.corners[i][3]; r2 = &mesh.corners[i][6]; l1 = &mesh.colors[i][3]; l2 = &mesh.colors[i][6]; break;
      }
     if(sqrt((r2[0]-r1[0])*(r2[0]-r1[0]) + (r2[1]-r1[1])*(r2[1]-r1[1]) + (r2[2]-r1[2])*(r2[2]-r1[2])) > 3E-6)
      {
      fprintf(outfile, "cylinder\n {\n");
      fprintf(outfile, " <%15f,%15f,%15f>, <%15f,%15f,%15f>,%15f\n", r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], rad);
      fprintf(outfile, " texture {\n  pigment { gradient <%f, %f, %f> scale %f translate <%f, %f, %f>\n",
              r2[0]-r1[0], r2[1]-r1[1], r2[2]-r1[2],
              sqrt((r2[0]-r1[0])*(r2[0]-r1[0]) + (r2[1]-r1[1])*(r2[1]-r1[1]) + (r2[2]-r1[2])*(r2[2]-r1[2])),
              r1[0], r1[1], r1[2]);
      if(transparency>1E-3) fprintf(outfile, "            color_map { [0, rgb <%.3f %.3f %.3f> transmit %.3f] [1, rgb <%.3f %.3f %.3f> transmit %.3f] } }",
                                    l1[0], l1[1], l1[2], transparency, l2[0], l2[1], l2[2], transparency);
       else fprintf(outfile, "            color_map { [0, rgb <%.3f %.3f %.3f>] [1, rgb <%.3f %.3f %.3f>] } }",
                    l1[0], l1[1], l1[2], l2[0], l2[1], l2[2]);
      if(finish==1) fprintf(outfile, "\n  finish { ambient BSAMBI diffuse BSDIFF specular BSSPEC } }\n }\n");
       else fprintf(outfile, " }\n }\n");
      }
     }
    }
   }
   else
    {
    for(i=0;i<mesh.ntriangles;i++)
     {
     for(j=0;j<3;j++)
      {
      switch(j)
       {
       case 0: r1 = &mesh.corners[i][0]; r2 = &mesh.corners[i][3]; break;
       case 1: r1 = &mesh.corners[i][0]; r2 = &mesh.corners[i][6]; break;
       case 2: r1 = &mesh.corners[i][3]; r2 = &mesh.corners[i][6]; break;
       }
      if(sqrt((r2[0]-r1[0])*(r2[0]-r1[0]) + (r2[1]-r1[1])*(r2[1]-r1[1]) + (r2[2]-r1[2])*(r2[2]-r1[2])) > 3E-6)
       {
       fprintf(outfile, "cylinder\n {\n");
       fprintf(outfile, " <%15f,%15f,%15f>, <%15f,%15f,%15f>,%15f\n", r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], rad);
       fprintf(outfile, " texture { pigment { rgb <%7.3f%7.3f%7.3f>", rgb[0], rgb[1], rgb[2]);
       if(transparency>1E-3) fprintf(outfile, " transmit %.3f }", transparency);
        else fprintf(outfile, " }");
       if(finish==1) fprintf(outfile, " finish { ambient BSAMBI diffuse BSDIFF specular BSSPEC } }\n }\n");
        else fprintf(outfile, " }\n }\n");
       }
      }
     }
    }
  fprintf(outfile, "\n");
  break;
 case 2:
  fprintf(outfile, "# mesh\nShape\n {\n");
  fprintf(outfile, " appearance Appearance { material charmolMaterial {");
  if(mesh.colors==NULL) fprintf(outfile, " emissiveColor %.3f %.3f %.3f", rgb[0], rgb[1], rgb[2]);
  if(transparency>1E-3) fprintf(outfile, " transparency %.3f } }\n", transparency);
   else fprintf(outfile, " } }\n");
  fprintf(outfile, " geometry IndexedLineSet\n  {\n");
  fprintf(outfile, "  coord Coordinate\n   {\n   point\n    [\n");
  for(i=0;i<mesh.ntriangles;i++)
   {
   fprintf(outfile, "%16f%16f%16f,\n%16f%16f%16f,\n%16f%16f%16f,\n",
           mesh.corners[i][0], mesh.corners[i][1], mesh.corners[i][2],
           mesh.corners[i][3], mesh.corners[i][4], mesh.corners[i][5],
           mesh.corners[i][6], mesh.corners[i][7], mesh.corners[i][8]);
   }
  fprintf(outfile, "    ]\n   }\n");
  fprintf(outfile, "  coordIndex\n   [\n");
  for(i=0,j=0;i<mesh.ntriangles;i++)
   {
   fprintf(outfile, "%16d%16d%16d%16d,\n", j, j+1, j+2, -1);
   j += 3;
   }
  fprintf(outfile, "   ]\n");
  if(mesh.colors!=NULL)
   {
   fprintf(outfile, "  color Color\n   {\n   color\n    [\n");
   for(i=0;i<mesh.ntriangles;i++)
    {
    fprintf(outfile, "%16f%16f%16f,\n%16f%16f%16f,\n%16f%16f%16f,\n",
            mesh.colors[i][0], mesh.colors[i][1], mesh.colors[i][2],
            mesh.colors[i][3], mesh.colors[i][4], mesh.colors[i][5],
            mesh.colors[i][6], mesh.colors[i][7], mesh.colors[i][8]);
    }
   fprintf(outfile, "    ]\n   }\n");
   }
  fprintf(outfile, "  }\n }\n\n");
  break;
 }

return outfile;
}

