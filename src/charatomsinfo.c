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

#include "charatomsinfo.h"

char **atomsinfo_fill_default_settings()
{
int rad, nrad, nsl;
char **atoms;

atoms = NULL;
nrad = 120;
nsl = 100;

if((atoms = (char **) malloc(nrad*sizeof(char *)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);}
rad = 0;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " X     0.30     0.3  0.3  0.3      0.0     0"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " H     0.31     1.0  1.0  1.0      1.0     1"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "He     0.28     0.5  0.5  0.5      4.0     2"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Li     1.28     0.5  0.5  0.5      6.9     3"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Be     0.96     0.5  0.5  0.5      9.0     4"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " B     0.84     0.0  0.8  0.8     10.8     5"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " C     0.76     0.7  0.2  0.0     12.0     6"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " N     0.71     0.0  0.8  1.0     14.0     7"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " O     0.66     1.0  0.0  0.0     16.0     8"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " F     0.57     0.5  1.0  0.2     19.0     9"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ne     0.58     0.5  0.5  0.5     20.2    10"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Na     1.66     0.5  0.5  0.5     23.0    11"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Mg     1.41     0.5  0.5  0.5     24.3    12"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Al     1.21     0.0  0.8  0.8     27.0    13"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Si     1.11     0.4  0.4  0.4     28.1    14"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " P     1.07     0.7  0.1  0.4     31.0    15"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " S     1.05     1.0  1.0  0.0     32.1    16"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Cl     1.02     0.0  0.8  0.0     35.5    17"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ar     1.06     0.5  0.5  0.5     39.9    18"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " K     2.03     0.5  0.5  0.5     39.1    19"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ca     1.76     0.5  0.5  0.5     40.1    20"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Sc     1.70     0.5  0.5  0.5     45.0    21"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ti     1.60     0.5  0.5  0.5     47.9    22"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " V     1.53     0.5  0.5  0.5     50.9    23"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Cr     1.39     0.5  0.5  0.5     52.0    24"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Mn     1.50     0.5  0.5  0.5     54.9    25"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Fe     1.42     0.5  0.5  0.5     55.8    26"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Co     1.46     0.5  0.5  0.5     58.9    27"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ni     1.24     0.5  0.5  0.5     58.7    28"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Cu     1.32     0.5  0.5  0.5     63.5    29"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Zn     1.22     0.5  0.5  0.5     65.4    30"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ga     1.22     0.0  0.8  0.8     69.7    31"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ge     1.20     0.7  0.7  0.7     72.6    32"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "As     1.19     0.7  0.1  0.4     74.9    33"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Se     1.20     1.0  0.6  0.0     79.0    34"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Br     1.20     0.7  0.1  0.1     79.9    35"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Kr     1.16     0.5  0.5  0.5     83.8    36"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Rb     2.20     0.5  0.5  0.5     85.5    37"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Sr     1.95     0.5  0.5  0.5     87.6    38"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " Y     1.90     0.5  0.5  0.5     88.9    39"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Zr     1.75     0.5  0.5  0.5     91.2    40"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Nb     1.64     0.5  0.5  0.5     92.9    41"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Mo     1.54     0.5  0.5  0.5     95.9    42"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Tc     1.47     0.5  0.5  0.5     98.0    43"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ru     1.46     0.5  0.5  0.5    101.1    44"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Rh     1.42     0.5  0.5  0.5    102.9    45"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Pd     1.39     0.5  0.5  0.5    106.4    46"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ag     1.45     0.5  0.5  0.5    107.9    47"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Cd     1.44     0.5  0.5  0.5    112.4    48"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "In     1.42     0.0  0.8  0.8    114.8    49"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Sn     1.39     0.7  0.7  0.7    118.7    50"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Sb     1.39     0.7  0.1  0.4    121.8    51"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Te     1.38     1.0  0.6  0.0    127.6    52"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " I     1.39     0.8  0.0  0.8    126.9    53"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Xe     1.40     0.5  0.5  0.5    131.3    54"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Cs     2.44     0.5  0.5  0.5    132.9    55"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ba     2.15     0.5  0.5  0.5    137.3    56"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "La     2.07     0.5  0.5  0.5    138.9    57"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ce     2.04     0.5  0.5  0.5    140.1    58"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Pr     2.03     0.5  0.5  0.5    140.9    59"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Nd     2.01     0.5  0.5  0.5    144.2    60"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Pm     1.99     0.5  0.5  0.5    145.0    61"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Sm     1.98     0.5  0.5  0.5    150.4    62"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Eu     1.98     0.5  0.5  0.5    152.0    63"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Gd     1.96     0.5  0.5  0.5    157.2    64"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Tb     1.94     0.5  0.5  0.5    158.9    65"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Dy     1.92     0.5  0.5  0.5    162.5    66"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ho     1.92     0.5  0.5  0.5    164.9    67"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Er     1.89     0.5  0.5  0.5    167.3    68"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Tm     1.90     0.5  0.5  0.5    168.9    69"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Yb     1.87     0.5  0.5  0.5    173.0    70"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Lu     1.87     0.5  0.5  0.5    175.0    71"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Hf     1.75     0.5  0.5  0.5    178.5    72"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ta     1.70     0.5  0.5  0.5    180.9    73"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " W     1.62     0.5  0.5  0.5    183.8    74"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Re     1.51     0.5  0.5  0.5    186.2    75"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Os     1.44     0.5  0.5  0.5    190.2    76"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ir     1.41     0.5  0.5  0.5    192.2    77"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Pt     1.36     0.5  0.5  0.5    195.1    78"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Au     1.36     0.5  0.5  0.5    197.0    79"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Hg     1.32     0.5  0.5  0.5    200.6    80"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Tl     1.45     0.0  0.8  0.8    204.4    81"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Pb     1.46     0.7  0.7  0.7    207.2    82"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Bi     1.48     0.7  0.1  0.4    209.0    83"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Po     1.40     1.0  0.6  0.0    210.0    84"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "At     1.50     0.5  0.5  0.5    210.0    85"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Rn     1.50     0.5  0.5  0.5    222.0    86"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Fr     2.60     0.5  0.5  0.5    223.0    87"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ra     2.21     0.5  0.5  0.5    226.0    88"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Ac     2.15     0.5  0.5  0.5    227.0    89"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Th     2.06     0.5  0.5  0.5    232.0    90"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Pa     2.00     0.5  0.5  0.5    231.0    91"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], " U     1.96     0.5  0.5  0.5    238.0    92"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Np     1.90     0.5  0.5  0.5    237.0    93"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Pu     1.87     0.5  0.5  0.5    244.0    94"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Am     1.80     0.5  0.5  0.5    243.0    95"); rad++;
if((atoms[rad] = (char *) malloc(nsl*sizeof(char)))==NULL) {puts("\nMemory allocation failed!\n"); exit(1);} strcpy(atoms[rad], "Cm     1.69     0.5  0.5  0.5    247.0    96"); rad++;
atoms[rad] = NULL;

return atoms;
}

