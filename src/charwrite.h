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

void write_run_system_command(char *inputcommand);

void write_determine_vrml_rotation(double *start, double *end, double *rot);

FILE *write_sphere(double *center, double rad, double *rgb, int finish, FILE *outfile, int outformat);

FILE *write_cylinder(double *start, double *end, double *rot, double rad, double *rgb, int finish, FILE *outfile, int outformat);

FILE *write_cylinder_topbot(double *start, double *end, double *rot, double rad, double *rgb, int finish, FILE *outfile, int outformat);

FILE *write_cone(double *start, double *end, double *rot, double rad, double *rgb, int finish, FILE *outfile, int outformat);

FILE *write_basic_settings(ATOMS *xyz, int natoms, TYCKA *tyc, KONUS *kon, int narrows, int camspec, int useangle, double *camera,
                           double convertxyzval, int scenerot, double scenerotangle, int camrot, double *camerarot, int camrq,
                           double *camerarq, int camsky, double *camerasky, int camzoom, double camerazoom, int handedness,
                           int finish, double *finset, FILE *outfile, int outformat);

FILE *write_xyz(ATOMS *xyz, int natoms, int nmolecules, int hideatoms, int atomcolor, int bondcolor, int finish, double *atomradiusscale,
                 double *lightatomradiusscale, double *bondradius, double **defatomcolor, double **defbondcolor, FILE *outfile, int outformat);

FILE *write_gauges(ATOMS *xyz, int natoms, GAUGE *gauge, int ngauges, int hideatoms, int finish, double *atomradiusscale,
                    double *lightatomradiusscale, double *bondradius, FILE *outfile, int outformat);

FILE *write_single_arrow(TYCKA tyc, KONUS kon, int finish, FILE *outfile, int outformat);

FILE *write_arrows(TYCKA *tyc, KONUS *kon, int narrows, int finish, FILE *outfile, int outformat);

FILE *write_mesh_solid(MESH mesh, double *rgb, double transparency, double isoval, int finish, FILE *outfile, int outformat);

FILE *write_mesh_wireframe(MESH mesh, double *rgb, double transparency, double isoval, int finish, FILE *outfile, int outformat);

