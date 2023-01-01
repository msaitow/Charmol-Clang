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

void check_basic_settings(int camspec, double *camera, int scenerot, double scenerotangle, int camrot, double *camerarot, int camrq,
                          double *camerarq, int camsky, double *camerasky, int camzoom, double camerazoom, int finish, double *finset);

void check_xyz(ATOMS *xyz, int natoms, int nmolecules, int hideatoms, int atomcolor, int bondcolor, double *atomradiusscale,
               double *lightatomradiusscale, double *bondradius, double **defatomcolor, double **defbondcolor);

void check_gauges(ATOMS *xyz, int natoms, GAUGE *gauge, int ngauges);

void check_single_arrow(TYCKA tyc, KONUS kon, int arroword, int vib);

void check_mesh(double *rgb, double transparency, double isoval, int surf);

