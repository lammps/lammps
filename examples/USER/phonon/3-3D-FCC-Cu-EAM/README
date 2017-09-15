This directory illustrates the usage of fix-phonon to calculate the dynamical
matrix as well as phonon dispersion curve for FCC Cu based on EAM potential.

The files under this directory:

 1) CuPhonon.bin.6500000   : last output binary file by fix-phonon
 2) CuPhonon.log           : log file for fix-phonon
 3) cuu3.eam               : EAM potential file for Cu
 4) data.pos               : LAMMPS input file
 5) disp.dat               : phonon dispersion data from CuPhonon.bin.6500000
 6) disp-expr.dat          : experimental phonon dispersion data for Cu
 7) disp-ld.dat            : phonon dispersion data by lattice dynamics based on EAM
 8) dos.dat                : phonon DOS data from CuPhonon.bin.6500000
 9) dos-expr.dat           : experimental PDOS for Cu
10) dos-ld.dat             : PDOS by LD based on EAM
11) in.disp/in.disp2       : input file to get disp.dat by phana
12) in.dos                 : input file to get dos.dat by phana
13) in.EAM3D               : LAMMPS input file
14) log.lammps             : LAMMPS log file
15) map.in                 : LAMMPS input file for fix-phonon
16) pdisp.eps              : figure of phonon dispersion curves
17) pdos.eps               : figure of phonon density of states
18) plot.disp              : gnuplot script to generate pdisp.eps
19) plot.dos               : gnuplot script to generate pdos.eps
20) pdisp.gnuplot          : gnuplot script to generate pdisp.eps (auto generated)
21) README                 : this file

To run this example, simply invoke: 
-> lmp -in in.EAM3D -screen none

Once done, one can use the auxiliary analysing code "phana" to obtain "disp.dat" and
"dos.dat" based on data from CuPhonon.bin.6500000:
-> phana CuPhonon.bin.6500000 < in.disp
-> phana CuPhonon.bin.6500000 < in.dos

And then use the gnuplot script file "plot.disp"/"plot.dos" to generate pdisp.eps/pdos.eps:
-> gnuplot plot.pdisp
-> gnuplot plot.pdos

The resultant ``pdisp.eps/pdos.eps'' compares the measured phonon dispersion to
experimental data and those by traditional lattice dynamics.

Alternatively, one can also use:
-> phana CuPhonon.bin.6500000 < in.disp2
-> gnuplot pdisp.gnuplot
to generate the phonon dispersion automatically.

NOTE: the binary file provided here might be unreadable on some computers because of
      incompatibility between different architectures.

Author: Ling-Ti Kong, konglt@sjtu.edu.cn
Nov 2015
