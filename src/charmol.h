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

static const double angstrom_to_bohr = 1.889725989;

typedef struct {char symb[20]; int atomnum; double coor[3]; int *cone; double rgb[3]; double rad; double scale; double mass; int molid; int ordid;} ATOMS;

typedef struct {int a1; int a2; double scale; double rad; int style; double rgb[3];} GAUGE;

typedef struct {char symb[20]; double scale; double start[3]; double end[3]; double rgb[3]; double rad; int ordid;} TYCKA;

typedef struct {double end[3]; double rad;} KONUS;

typedef struct {int number; double scale; int atomiccolor; double rgb[3]; double arradscale; double conradscale; int *atoms;} VIBRATION;

typedef struct {int atom; char type[5]; double coor[3]; int degen; int nprim; double *expon; double *contr; double *contr2; double screenr2;} BASIS;

typedef struct {char label[20]; int spin; double *coef;} MOS;

typedef struct {int order; char label[20]; int spin; double isoval; int gencon; int screen; double screenfactor; double screenval; double gridsizescale;
                char griddensity[20]; double griddist[3]; double transparency; double rgbpos[3]; double rgbneg[3]; double *coef; int longname; int style;} ORBITAL;

typedef struct {char filename[500]; FILE *file; int norb; int *orbs;} ORBCOMB;

typedef struct {char filename[500]; double isoval; int sign; double transparency; double rgb[3]; double rgbpos[3]; double rgbneg[3];
                int nsurf; int *surfs; int norb; int *orbs; int longname; int style; int colmap; char potfilename[500];} SURFACE;

typedef struct {int ntriangles; double **corners; double **normals; double **colors;} MESH;

