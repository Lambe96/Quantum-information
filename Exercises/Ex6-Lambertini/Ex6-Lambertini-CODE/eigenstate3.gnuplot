set datafile separator ','
set grid xtics mxtics ytics mytics back
set title '𝜓_1(x) comparison'
set xlabel 'x'
set ylabel '𝜓(x)'
set key right bottom
set xrange [-5:5]
f(x) = (1/(2*pi))**(0.25)*x*exp(-x**2/4)

plot 'x_and_columns.txt' u 1:3 w l linewidth 3 lc 'red' title 'computed 𝜓_1(x) ', f(x) w l linewidth 3 dashtype 3 lc 2 title 'theoretical 𝜓_1(x)'
