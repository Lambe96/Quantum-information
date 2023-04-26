set datafile separator ','
set grid xtics mxtics ytics mytics back
set title 'Groundstates for different values of  ƛ computed with RSRG'
set xlabel 'ƛ'
set ylabel 'E_0(ƛ)'
set key left bottom
set xrange [-0.2:3.2]
set yrange [-3.5:-0.8]

f(x) = -1-x**2/4
g(x) = -x
h(x) = (x<2)?f(x):g(x)


 plot   'ground_states.txt'  w l title sprintf("E_0"),[0:3] h(x) lc 'black' title 'Mean field'
