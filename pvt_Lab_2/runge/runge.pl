#!/usr/bin/gnuplot
#!/usr/bin/gnuplot -persist

set terminal png size 1200,800 linewidth 4 font "Verdana,14"
set termoption dash
set output "runge.png"
set grid x y
set xlabel "Number of processes"
set ylabel "Relative speedup (2 processes)"
set key left top
unset log
set xrange [0 : 70]
set yrange [0 : 70]
plot "linear.txt" using 1:2 with linespoints linecolor 4 title "Linear speedup", "out.txt" using 1:2 with linespoints linecolor 5 title "Parallel speedup"
