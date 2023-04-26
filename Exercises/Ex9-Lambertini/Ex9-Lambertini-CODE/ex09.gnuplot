set datafile separator ','
set grid xtics mxtics ytics mytics back
set title 'Groundstates for systems with N spins'
set xlabel 'ƛ'
set ylabel 'E_0(ƛ)'
set key left bottom
set xrange [-0.2:3.2]
set yrange [-7:-0.8]

f(x) = -1-x**2/4
g(x) = -x
h(x) = (x<2)?f(x):g(x)


 plot for [k = 2:12]  'gs_N.txt' u 1:k w l title sprintf("E_0, N=".(k)),[0:3] h(x) lc 'black' title 'Mean field'
