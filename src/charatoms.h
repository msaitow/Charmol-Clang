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

double atoms_check_xyz_units(int extfileformat, char **input, char **extinput);

void atoms_rescale_xyz(ATOMS *xyz, int natoms, double scaleval);

void atoms_fill_basic_settings(ATOMS *xyz, int natoms, char **atoms);

void atoms_change_atomscales(ATOMS *xyz, int natoms, char **input);

void atoms_change_atomcolors(ATOMS *xyz, int natoms, char **input);

void atoms_make_connectivitymatrix(ATOMS *xyz, int natoms, char **input);

void atoms_make_gauge_settings(GAUGE *gauge, int ngauges, ATOMS *xyz, int natoms, int nmolecules, double *bondradius, char **input);

