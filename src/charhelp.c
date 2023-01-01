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

#include "charhelp.h"

void help_print()
{

puts("");
puts("");
puts("                                                     CHARMOL");
puts("  ");
puts("                        program for making high-quality pictures of molecular structures");
puts("                             using POV-Ray (Persistence of Vision Raytracer) or VRML");
puts("");
puts("");
puts("syntax: charmol inputfile [povrayverbose] [hiderendering/norendering]");
puts("        - 'povrayverbose' option forces printing of information provided by povray during rendering");
puts("        - 'hiderendering' option supresses showing progress of povray rendering in graphical form");
puts("        - 'norendering' option supresses povray rendering");
puts("        - '-h' or '--help' options serve for printing the help");
puts("");
puts("1) charmol is a program for making high-quality pictures of molecular structures");
puts("   - it is capable of doing pictures of more molecules together, depiction of orbitals and other surfaces, schematic depiction of vibrations");
puts("     and other user defined arrows, measuring gauges, fine tuning of the design according to ones needs (user defined colors, radii), etc.");
puts("   - default settings are carefully chosen in order to produce reasonable results without need of 'advanced' input editing");
puts("   - simple ascii file serves as input, its output can be either vrml or povray file containing the 3D scene");
puts("     if povray output is chosen, povray is automatically called by the charmol program to create the picture");
puts("");
puts("2) structure of the input file");
puts("   2a) lines beginning with '#' are comments");
puts("   2b) input is normaly composed of sections and keywords with their arguments");
puts("   2c) every global (not belonging to a section) keyword must be on a separate line and outside all the sections (with 'cubefile' exception),");
puts("       otherwise some of the required user-specified settings could be neglected (possibly without any error message)");
puts("       keywords belonging to sections can be (has to be) arbitrarily combined on single line inside of the particular section");
puts("       reading of the input is case sensitive - all sections, keywords and string arguments has to be specified like they are printed here");
puts("       (i.e. charmol input file should be composed only of small letters except names of external files and possibly element symbols)");
puts("   2d) sections:");
puts("       - xyz [au/angstrom] (endxyz)     - mandatory, if none of the molden, fchk or cube external files was used");
puts("                                          if geometry is in atomic units, 'au' has to be specified on the line containing 'xyz' directive");
puts("                                          (by default, geometry is expected to be in angstroms)");
puts("                                          lines between xyz and endxyz are interpreted, each line must contain element symbol and xyz coordinates");
puts("                                          it can optionaly also contain molecular identification (mid) if you want more molecules in 1 picture");
puts("                                          atoms with the same value of user-specified integer mid are interpreted as one molecule");
puts("                                          when mids are used, make sure that some mid value is ascribed to each atom");
puts("                                          mids can be also specified in external molden file, cube file or fchk file");
puts("                                          (for mids in fchk, add section 'Molecular ids' with syntax like 'Atomic numbers' section has)");
puts("                                          syntax of a line: element_symbol x y z [mid]");
puts("                                          (examples: H 0.0 0.0 0.0, Cl 1.0 0.0 0.3 1)");
puts("       - atomcolor (endatomcolor)       - optional, lines between atomcolor and endatomcolor are interpreted, each line contains information");
puts("                                          about color (RGB in <0.0,1.0>) of atoms specified by either element symbol or number of atom");
puts("                                          (examples of a line: 2 0.2 0.3 0.2, Cl 0.0 1.0 0.0)");
puts("       - atomscale (endatomscale)       - optional, lines between atomscale and endatomscale are interpreted, each line contains information");
puts("                                          about scale of depicted radii of atoms specified by either element symbol or number of atom");
puts("                                          (examples of a line: 2 0.9, Cl 1.3)");
puts("       - arrow (endarrow)               - optional, lines between arrow and endarrow are interpreted, each line contains information");
puts("                                          about one arrow which should be added to the picture");
puts("                                          each line can have three different forms:");
puts("                                           1. atom_number scale x y z r g b");
puts("                                              - places arrow defined by (x,y,z) and its scaling factor 'scale' on atom 'atom_number'");
puts("                                                arrow has (r,g,b) color in RGB scale (each number in <0.0,1.0>)");
puts("                                           2. orig scale x y z r g b");
puts("                                              - places arrow to the origin of coordinate system");
puts("                                           3. coor scale x y z r g b i j k");
puts("                                              - places arrow to user specified coordinates (i,j,k)");
puts("                                          (example of a line: 9 0.5 0.000000 -0.257785 -3.818825 0.8 0.2 0.8)");
puts("       - vibration (endvibration)       - optional, lines between vibration and endvibration are interpreted");
puts("                                          depiction of vibrations works with external molden file");
puts("                                          lines not containing vibration specification change settings for all vibrations,");
puts("                                          otherwise settings are changed for particular vibration (overwriting previous settings)");
puts("                                          each line can contain following keywords and their arguments:");
puts("                                           1. vib n                     - specifies the vibration to depict (n is its order number: 1-3N)");
puts("                                           2. vibrationscale x          - changes length of arrows (default x=1.0)");
puts("                                           3. arradscale x              - changes radius of arrows (default x=1.0)");
puts("                                                                          (conradscale is scaled accordingly)");
puts("                                           4. conradscale x             - changes radius of cones (default x=1.0)");
puts("                                           5. atomiccolors on/off       - turns on/off atomic colors for arrows (default off)");
puts("                                           6. vibrationcolor r g b      - changes RGB color of arrows (default black - r=0.0 g=0.0 b=0.0)");
puts("                                           7. includeatoms a1 a2 ...    - specifies atoms, for which displacement will be shown");
puts("                                           8. excludeatoms a1 a2 ...    - specifies atoms, for which displacement will not be shown");
puts("                                          (atoms in includeatoms and excludeatoms can be specified separately or by ranges aN-aM)");
puts("                                          (by default, displacements will be shown for all atoms of the molecule)");
puts("                                          (example of a line: vib 7 vibrationscale 0.8 excludeatoms 1 3-6 9)");
puts("       - orbital (endorbital)           - optional, lines between orbital and endorbital are interpreted");
puts("                                          depiction of orbitals works with external molden or fchk file and GTO basis sets");
puts("                                          lines not containing orbital specification change settings for all orbitals,");
puts("                                          otherwise settings are changed for particular orbitals (overwriting previous settings)");
puts("                                          each line can contain following keywords and their arguments:");
puts("                                           1. generalcontraction on/off - turns on/off using general contraction");
puts("                                                                          (default on for generally-contracted basis sets)");
puts("                                           2. prescreening on/off       - turns on/off using prescreening of contributions of GTOs (default on)");
puts("                                                                          if you are not sure about the result, turn the prescreening off");
puts("                                           3. prescreeningfactor x      - changes factor used for prescreening of contributions of GTOs");
puts("                                                                          contributions smaller than 'x*isovalue' are neglected (default x=0.001)");
puts("                                           4. gridsizescale x           - rescales the grid size used for calculation (default x=1.0)");
puts("                                                                          (useful especially for very diffuse orbitals, for example Rydbergs)");
puts("                                           5. griddensity density       - changes the density of the grid");
puts("                                                                          density can be: verylow, low, medium, high, veryhigh (default medium)");
puts("                                           6. longfilenames on/off      - turns on/off using long filenames (default off)");
puts("                                                                          can be useful for tuning the settings for depiction of orbitals");
puts("                                           7. isovalue x                - changes isovalue used for depiction (default x=0.05)");
puts("                                           8. positivecolor r g b       - changes color of positive part of orbital (default blue - r=0.0 g=0.0 b=1.0)");
puts("                                           9. negativecolor r g b       - changes color of negative part of orbital (default red - r=1.0 g=0.0 b=0.0)");
puts("                                          10. transparency x            - changes transparency of orbital (default no transparency - x=0.0)");
puts("                                          11. style type                - changes the style of orbital (type can be 'solid' or 'wireframe', default solid)");
puts("                                          12. alpha o1 o2 ...           - specifies alpha orbitals to calculate and depict");
puts("                                          13. beta o1 o2 ...            - specifies beta orbitals to calculate and depict");
puts("                                          14. combineorbitals alpha o1 o2 ... beta o1 o2 ...");
puts("                                                                        - produces an extra file containing combination of orbitals");
puts("                                                                          (orbitals settings are taken from previous specifications for particular orbitals)");
puts("                                          (orbitals can be specified separately or by ranges oN-oM, each oN can be specified in two different ways:");
puts("                                           if symmetry of orbital is used, specification is taken as an orbital label, otherwise it is taken as an");
puts("                                           orbital ordering number - for example '3a1' and '7' could represent the same orbital)");
puts("                                          (if you use 'all' as o1 argument, all alpha or beta orbitals will be done)");
puts("                                          (example of a line: alpha 1-2 beta 1 isovalue 0.05 transparency 0.3)");
puts("       - surface (endsurface)           - optional, lines between surface and endsurface are interpreted");
puts("                                          depiction of surfaces works with external gaussian cube files with orthogonal voxels and");
puts("                                          data for scalar properties like densities and orbitals");
puts("                                          there can be specified arbitrary number of cube files (multiple cubefile directives)");
puts("                                          geometry is read from the first cube file found (it can be either outside or inside the surface section)");
puts("                                          only cube files specified inside the surface section are used for depiction of surfaces");
puts("                                          lines not containing cubefile specification change settings for all surfaces,");
puts("                                          otherwise settings are changed for particular surfaces (overwriting previous settings)");
puts("                                          each line can contain following keywords and their arguments:");
puts("                                           1. cubefile filename         - specifies the gaussian cube file to use");
puts("                                           2. isovalue x                - changes isovalue used for depiction");
puts("                                                                          (default x=0.05 for orbitals and x=0.1 for other types of surfaces)");
puts("                                                                          (surface type is determined automatically from the cube file)");
puts("                                           3. longfilenames on/off      - turns on/off using long filenames (default off)");
puts("                                                                          can be useful for tuning the settings for depiction of surfaces");
puts("                                           4. color r g b               - changes color of non-orbital surface (default yellow - r=1.0 g=1.0 b=0.0)");
puts("                                           5. positivecolor r g b       - changes color of positive part of orbital (default blue - r=0.0 g=0.0 b=1.0)");
puts("                                           6. negativecolor r g b       - changes color of negative part of orbital (default red - r=1.0 g=0.0 b=0.0)");
puts("                                           7. transparency x            - changes transparency of surface (default no transparency - x=0.0)");
puts("                                           8. style type                - changes the style of surface (type can be 'solid' or 'wireframe', default solid)");
puts("                                           9. colormapped filename      - additional gaussian cube file (containing potential) is used to make color-mapped suface");
puts("                                          10. orbitals o1 o2 ...        - specifies orbitals to depict for cube file containing more orbitals");
puts("                                                                          (all orbitals found in cube file will be done by default)");
puts("                                          (orbitals can be specified separately or by ranges oN-oM)");
puts("                                          (example of a line: cubefile test.cube isovalue 0.5 transparency 0.2)");
puts("");
puts("   2e) keywords:");
puts("       - keywords connected to making the 3D scene:");
puts("         - addbond a1 a2                - adds bond connecting atoms a1 and a2 (atom numbers)");
puts("                                          for the case that bond does not appear in the picture");
puts("                                          (bonds are determined automatically by calculation of connectivity matrix from covalent radii)");
puts("         - delbond a1 a2                - deletes bond connecting atoms a1 and a2 (atom numbers)");
puts("         - delallbonds a1               - deletes all bonds connected to atom a1 (atom number)");
puts("         - addgauge a1 a2 [x s r g b]   - adds measuring gauge connecting atoms a1 and a2 (atom numbers)");
puts("                                          gauge is a special bond for measuring the distance between atoms");
puts("                                          optional 'x s r g b' has following meaning and defaults:");
puts("                                          x - scale of gauge radius (radius is 'x' times radius of a real bond, default x=0.4)");
puts("                                          s - style of gauge, it can be dashed, dotted, dotdashed or solid (default s=dashed)");
puts("                                          r g b - color of gauge (default black - r=0.0 g=0.0 b=0.0)");
puts("         - atomradiusscale x [mid]      - changes the scale of spheres representing atoms (default x=0.4)");
puts("                                          thickness of other objects (bonds and arrows) is scaled accordingly");
puts("                                          when optional mid is used, settings are changed only for corresponding molecule");
puts("                                          when mid is not used, it changes settings for all molecules present in the input");
puts("                                          (independently on the fact, whether mid was or was not used in the xyz section)");
puts("         - lightatomradiusscale x [mid] - changes the scale of spheres representing light atoms (default x=1.5)");
puts("                                          (see section xyz and keyword atomradiusscale for details about meaning of optional mid)");
puts("         - bondradiusscale x [mid]      - changes the scale of cylinders representing bonds (default x=1.0)");
puts("                                          (see section xyz and keyword atomradiusscale for details about meaning of optional mid)");
puts("         - arrowradiusscale x [mid]     - changes the scale of cylinders representing arrows (default x=1.0)");
puts("                                          thickness of cones on the end of arrows is scaled accordingly");
puts("                                          (see section xyz and keyword atomradiusscale for details about meaning of optional mid)");
puts("         - coneradiusscale x [mid]      - changes the scale of cones on the end of arrows (default x=1.0)");
puts("                                          (see section xyz and keyword atomradiusscale for details about meaning of optional mid)");
puts("         - hideatoms                    - supresses the atomic spheres depiction (no argument)");
puts("         - noatomcolor                  - supresses using different colors for different elements (no argument)");
puts("                                          when specified defaultatomcolor is automatically used");
puts("         - nobondcolor                  - supresses the bicolor bonds depiction (no argument)");
puts("                                          when specified defaultbondcolor is automatically used");
puts("         - defaultatomcolor r g b [mid] - changes the default atom color (default black - r=0.0 g=0.0 b=0.0)");
puts("                                          when specified noatomcolor is automatically used");
puts("                                          if neither nobondcolor nor defaultbondcolor is used, defaultatomcolor is used also for bonds");
puts("                                          (see section xyz and keyword atomradiusscale for details about meaning of optional mid)");
puts("         - defaultbondcolor r g b [mid] - changes the default bond color (default grey - r=0.7 g=0.7 b=0.7)");
puts("                                          when specified nobondcolor is automatically used");
puts("                                          (see section xyz and keyword atomradiusscale for details about meaning of optional mid)");
puts("         - nocenterofmass               - supresses moving the center of mass of the molecule to the origin of coordinate system");
puts("                                          (no argument, when used depiction is made using original coordinates)");
puts("         - atomsinfofile file           - user specified file with atomic properties is used instead of default settings");
puts("                                          file should contain line per atom with syntax 'element_symbol covalent_radius r g b mass atomic_number'");
puts("                                          (i.e. symbol of the element, its covalent radius, color in RGB scale, its atomic mass and number of protons)");
puts("                                          (default settings can be found in {your_destination_of_release}/src/charatomsinfo.info)");
puts("         - moldenfile file              - user specified molden file is used, only one (first specified) molden file is used");
puts("                                          (atomic colors, scales, arrows, gauges etc. specified in charmol input are still used)");
puts("         - fchkfile file                - user specified fchk file is used, only one (first specified) fchk file is used");
puts("                                          (atomic colors, scales, arrows, gauges etc. specified in charmol input are still used)");
puts("         - cubefile file                - user specified gaussian cube file is used, only one (first specified) cube file is used for");
puts("                                          reading the geometry (it can be, but does not have to be, specified outside the surface section)");
puts("                                          see section surface for more details about using gaussian cube files");
puts("                                          (atomic colors, scales, arrows, gauges etc. specified in charmol input are still used)");
puts("         - output format                - specification of output file format, format can be povray or vrml (default is povray)");
puts("                                          in case of vrml, its version 2.0 is used (also known as VRML97)");
puts("");
puts("       - keywords connected to povray settings (see povray documentation for detailed explanation of some of them):");
puts("         (note that some of these keywords automatically supress using of some others to avoid meaningless transfomations)");
puts("         - cameraposition x y z         - changes the location of camera to (x,y,z)");
puts("                                          by default location of camera is automatically calculated along the z axis");
puts("         - camerarotquat x y z w        - rotates the camera position according to specified rotational quaternion");
puts("                                          meaning of the parameters come from the quaternion definition q = w + ix + jy + kz");
puts("         - camerazoom value             - performs zooming by value in percents, positive values correspond to zoom in");
puts("         - scenerotation angle          - rotates the scene around vector of view (picture in its own plane) by angle in degrees");
puts("                                          positive values correspond to clockwise rotation, independently on handedness chosen");
puts("         - camerasky x y z              - changes the camera sky vector (povray default is (0,1,0))");
puts("                                          this can serve for the same thing as scenerotation");
puts("         - camerarotation x y z         - rotates the camera around the 3 axes by (x,y,z) in degrees");
puts("                                          serves for specification of general rotation of the scene");
puts("         - handedness right/left        - changes the scene handedness (default is right-handed scene)");
puts("         - noangleuse                   - supresses usage of camera angle (default 15° angle is used)");
puts("         - moldenview                   - invokes production of scene from the same view like molden uses");
puts("                                          (scene is viewed from z axis, right handed coordinate system and scenerotation by -90° is used)");
puts("         - finish i j k                 - changes settings for finishing the design (default i=0.2 j=0.8 k=0.8)");
puts("                                          (i stands for ambient, j for diffuse and k for specular)");
puts("         - nofinish                     - supresses finishing of the design (finishing is used by default)");
puts("");
puts("       - keywords connected to vrml settings (see vrml documentation for detailed explanation):");
puts("         - ambient i                    - changes the ambient intensity (default i=0.2)");
puts("         - shininess i                  - changes the reflection shininess (default i=0.25)");
puts("         - nofinish                     - supresses finishing of the design (finishing is used by default)");
puts("");
puts("3) setting the POV-Ray camera");
puts("   - there are three ways how to set the povray camera with charmol program, they are listed here from the most to");
puts("     the least 'easy-to-use' one (3a and 3b may alter depending on the case):");
puts("   3a) use 'charcam' program (part of the charmol release)");
puts("       - charcam is an interactive script which serves for setting the povray camera");
puts("       - it is capable of making rotations around 3 axes (horizontal, vertical, and axis perpendicular to the plane");
puts("         of monitor) and zooming");
puts("       - all transformations done by the script are performed in the 'actual' coordinate system (like you would expect");
puts("         when rotating the scene with mouse) using general rotations described by rotational quaternions");
puts("       - try 'charcam -h' for more information");
puts("   3b) use camera settings provided by VRML player (or other program capable of working with vrml or povray in 3D):");
puts("       - camera settings can be represented either by rotational quaternion (use 'camerarotquat' charmol keyword),");
puts("         or by camera position and sky vectors (use 'cameraposition' and 'camerasky' keywords)");
puts("       - since vie3dscene is free and good VRML player for linux, we will show here the procedure with view3dscene in");
puts("         more details:");
puts("         - produce vrml output by charmol ('output vrml'), open it with view3dscene and find the desired view of the scene");
puts("         - use 4-number vector from low-left corner with label 'Rotation quat' for 'camerarotquat' keyword (older versions of");
puts("           view3dscene, in which rotational quaternion information is printed, for example 3.5.2) or");
puts("           use 'Console -> Print Current Camera (VRML 2.0)' from menu and use printed vectors 'position' for 'cameraposition'");
puts("           keyword and 'up' for 'camerasky' keyword, note that additional zooming using 'camerazoom' might be needed (newer");
puts("           versions of view3dscene, for example 3.11)");
puts("         - produce povray output with the new camera settings by charmol");
puts("   3c) set the camera by hand by changing first 'cameraposition' keyword parameters and then tune the desired view by using");
puts("       'scenerotation' keyword (or any other combination of possibilities provided in charmol program, see the section of help");
puts("        about keywords connected to povray settings)");
puts("");
puts("4) OpenMP parallelization");
puts("   - time-consuming parts of the code are parallelized using OpenMP, you can control the OpenMP functionality in a standart way,");
puts("     for example number of processors to use is controled by setting the OMP_NUM_THREADS environment variable");
puts("");
puts("5) running the charmol program");
puts("   - except classical run from command line on local machine, charmol can be also easily run on remote hosts (that can be convenient");
puts("     for example in case of calculation of large number of orbitals), in such a case it is rather recommended to use version");
puts("     with statically linked libraries to avoid version-related problems (see file {your_destination_of_release}/INSTALL for details)");
puts("");
puts("6) input examples");
puts("   - some examples of input files are located in the directory {your_destination_of_release}/examples");
puts("");
puts("");
puts("For bug reporting and suggestions about additional functions, please contact the author.");
puts("   Jakub Chalupsky       (email: chalupsky.jakub@gmail.com)");
puts("");
puts("");

}

