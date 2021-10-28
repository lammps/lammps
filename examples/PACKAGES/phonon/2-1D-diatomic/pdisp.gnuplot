set term post enha colo 20
set out "pdisp.eps"

set xlabel "q"
set ylabel "frequency (THz)"

set xrange [0:1]
set yrange [0:*]

set grid xtics
# {/Symbol G} will give you letter gamma in the label
set xtics ("" 0, "" 1)

unset key

plot "disp.dat" u 4:5 w l lt 1,\
"" u 4:6 w l lt 1