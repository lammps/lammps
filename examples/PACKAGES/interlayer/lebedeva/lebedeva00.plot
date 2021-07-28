# Author: Zbigniew Koziol, National Center for Nuclear Research, Poland
# Email: softquake@gmail.com

set term x11;
unset log
unset title
set size 1.0,1.0

set encoding iso_8859_1
#set term postscript eps enhanced color "Helvetica" 18;
#set output "lebedeva00.eps"

set zero 1e-018;

set xlabel "x,y [{\305}]" font "Helvetica,18";
set ylabel "U [eV/atom]" font "Helvetica,18";

set key font ",18"

set key right
set key top

set pointsize 1.2

set xrange [0:20]
set yrange [-0.002:0.001]
#set yrange [-0.01:0.01]
#set yrange [*:*]

plot \
	"LamppsResult.dat" u 2:5 t "Leb LAMMPS",\
	"PerlResult.dat" u 1:($4*0.001/2.) w l t "Leb Perl"

exit
