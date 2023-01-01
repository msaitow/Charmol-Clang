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

void mesh_calc_gradients(int *gp, double *gd, double *gval, int gi, int gj, int gk, double **cubecorgrad);

void mesh_calc_colors(int *gp, double *gpot, int gi, int gj, int gk, double **cubecorcol);

MESH mesh_process_cube(int *gp, double *gd, double *gx, double *gy, double *gz, double *gval, double *gpot, double isoval, int gi, int gj, int gk,
                       double *cubecorval, double **cubecorcoor, double **cubecorgrad, double **cubecorcol, MESH mesh);

MESH mesh_clean_mesh(MESH mesh);

MESH mesh_generate_mesh(int *gp, double *gd, double *gx, double *gy, double *gz, double *gval, double *gpot, double isoval);

