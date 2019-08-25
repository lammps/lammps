set term post enha colo 20
set out "pdisp.eps"

set xlabel "q"
set ylabel "frequency (THz)"

set xrange [0:1.68817]
set yrange [0:*]

set grid xtics
# {/Symbol G} will give you letter gamma in the label
set xtics ("" 0, "" 0.707107, "" 0.942809, "" 1.68817)

unset key

plot "pdisp.dat" u 4:5 w l lt 1,\
"" u 4:6 w l lt 1,\
"" u 4:7 w l lt 1,\
"" u 4:8 w l lt 1,\
"" u 4:9 w l lt 1,\
"" u 4:10 w l lt 1