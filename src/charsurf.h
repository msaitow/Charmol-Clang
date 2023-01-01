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

void surf_process_options(SURFACE *surf, int num, char *line, int global);

void surf_make_surfaces(ATOMS *xyz, int natoms, int nmolecules, int centerofmass, double *centerofmasscoor, double convertxyzval,
                        int hideatoms, int atomcolor, int bondcolor, double *atomradiusscale, double *lightatomradiusscale,
                        double *bondradius, double **defatomcolor, double **defbondcolor, TYCKA *tyc, KONUS *kon, int narrows,
                        GAUGE *gauge, int ngauges, int camspec, int useangle, double *camera, int scenerot, double scenerotangle,
                        int camrot, double *camerarot, int camrq, double *camerarq, int camsky, double *camerasky, int camzoom,
                        double camerazoom, int handedness, int finish, double *finset, char **input, char **extinput,
                        int extfileformat, char *outfilename, char *suffix, int outformat, char *command, char *comend);

