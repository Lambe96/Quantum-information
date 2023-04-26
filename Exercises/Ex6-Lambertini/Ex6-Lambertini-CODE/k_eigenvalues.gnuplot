set datafile separator ','
set grid xtics mxtics ytics mytics back
set title 'Eigenvaues comparison'
set xlabel 'n'
set ylabel 'E(n)'
set key right bottom
set xrange [0:110]
f(x) = x+ 0.5

plot 'first_k_eigenvalues.txt' w l linewidth 3 lc 'red' title 'computed eigenvalues', f(x) w l linewidth 3 dashtype 3 lc 2 title 'theoretical eigenvalues'
