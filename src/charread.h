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

char **read_file_to_char(char *filename);

char **read_xyzpartofcubefile_to_char(char *filename);

char **read_file_to_char_withoutcomments(char *filename);

char **read_repair_fortran_single_precision(char **input, char *filename);

char *read_next_token(char *string);

char *read_firststrcmp(char *str, char *tok);

ATOMS *read_xyz_input(char **input);

ATOMS *read_xyz_molden(char **input);

ATOMS *read_xyz_fchk(char **input);

ATOMS *read_xyz_cube(char **input);

BASIS *read_gto_molden(ATOMS *xyz, int natoms, char **input);

BASIS *read_gto_fchk(ATOMS *xyz, int natoms, char **input);

MOS *read_mo_molden(int nbasis, char **input);

MOS *read_mo_fchk(int nbasis, char **input);

void read_check_geom_surf_cube(SURFACE *surf, int s, ATOMS *xyz, int natoms, int centerofmass, double *centerofmasscoor, char **input);

void read_grid_info_cube(int *gpoints, double *gmin, double *griddist, int centerofmass, double *centerofmasscoor, char **input);

void read_check_geom_grid_cube(ATOMS *xyz, int natoms, int centerofmass, double *centerofmasscoor, int *gpoints, double *gmin, double *griddist, char **input);

void read_grid_values_oneset_cube(double *gridval, int totpoints, char **input);

void read_grid_values_cube(double **gridval, int totpoints, int totsurfs, char **input);

void read_check_grid_value_signs_cube(SURFACE *surf, int s, double **gridval, int totpoints, int totsurfs);

void read_vibration_molden(TYCKA *vibtyc, KONUS *vibkon, ATOMS *xyz, int natoms, int centerofmass, double *centerofmasscoor, int vibnum, char **input);

