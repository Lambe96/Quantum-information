set datafile separator ','
set grid xtics mxtics ytics mytics back
set title 'First four eigenvalues for systems with 10 spins'
set xlabel 'ƛ'
set ylabel 'E_0(ƛ)'
set key left bottom




plot for [k = 2:5]   'firstfoureigv_10.txt' u 1:k w l title sprintf("E_N N=".(k))
