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

double calc_prim(double r2, int nprim, double *expon, double *contr);

double calc_prim_sum_gridexp(int nprim, double *contr, double *gridexp, int ge);

double calc_screen_factor_shell(int degen, int nprim, double *expon, double *contr, double *contr2, double r2);

void calc_screenr2_values(BASIS *gto, int natoms, int nshell, int maxdegen, double screenval, int gencon);

double calc_shell_molden(double r2, double x, double y, double z, int degen, int nprim, double *expon, double *contr, double *contr2,
                         double *motocalc, int mo, double *gridexp, int ge, int gencon);

double calc_shell_fchk(double r2, double x, double y, double z, int degen, int nprim, double *expon, double *contr, double *contr2,
                       double *motocalc, int mo, double *gridexp, int ge, int gencon);

void calc_grid_gencon_molden(BASIS *gto, int natoms, int nshell, int maxdegen, int totpoints, int screen, double *gridx, double *gridy, double *gridz, double *gridval,
                             double *motocalc, double *gridexp, int gencon);

void calc_grid_gencon_fchk(BASIS *gto, int natoms, int nshell, int maxdegen, int totpoints, int screen, double *gridx, double *gridy, double *gridz, double *gridval,
                           double *motocalc, double *gridexp, int gencon);

void calc_grid_segcon_molden(BASIS *gto, int natoms, int nshell, int maxdegen, int totpoints, int screen, double *gridx, double *gridy, double *gridz, double *gridval,
                             double *motocalc, double *gridexp, int gencon);

void calc_grid_segcon_fchk(BASIS *gto, int natoms, int nshell, int maxdegen, int totpoints, int screen, double *gridx, double *gridy, double *gridz, double *gridval,
                           double *motocalc, double *gridexp, int gencon);

int calc_check_gencon(int natoms, int nshell, int maxdegen, BASIS *gto);

