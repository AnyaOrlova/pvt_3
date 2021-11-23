#!/usr/bin/gnuplot
#!/usr/bin/gnuplot -persist
#!/bin/bash

set terminal png size 1200,800 linewidth 4 font "Verdana,14"
set termoption dash
set output "sgemv.png"
set grid x y
set xlabel "Number of processes"
set ylabel "Relative speedup (2 processes)"
set key right bottom
unset log
set xrange [8 : 64]
set xtics nomirror 8
set yrange [0 : 64]
plot "linear.txt" using 1:2 with linespoints linecolor 1 title "Linear speedup", "28000.txt" using 1:2 with linespoints linecolor 3 title "SGEMV speedup for n=m=28000", "45000.txt" using 1:2 with linespoints linecolor 2 title "SGEMV speedup for n=m=45000"
