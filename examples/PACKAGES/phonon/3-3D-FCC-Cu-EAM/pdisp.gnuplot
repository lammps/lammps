set term post enha colo 20
set out "pdisp.eps"

set xlabel "q"
set ylabel "frequency (THz)"

set xrange [0:4.47091]
set yrange [0:*]

set grid xtics
# {/Symbol G} will give you letter gamma in the label
set xtics ("{/Symbol G}" 0, "X" 0.707107, "W" 1.06066, "K" 1.23744, "{/Symbol G}" 2.156, "L" 3.02202, "U" 3.32821, "W" 3.50498, "L" 3.85854, "K/U" 4.16472, "X" 4.47091)

unset key

plot "disp.dat" u 4:5 w l lt 1,\
"" u 4:6 w l lt 1,\
"" u 4:7 w l lt 1