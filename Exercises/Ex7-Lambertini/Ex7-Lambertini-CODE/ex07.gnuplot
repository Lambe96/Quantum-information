set datafile separator ','
set grid xtics mxtics ytics mytics back
set title 'time evolution of the system from t=0 to t=3 with 1000 time frame'
set xlabel 'x'
set ylabel 'P(x)'
set key right bottom

 plot for [k = 2:1002:200]  'time_evolution.txt' u 1:k w l title sprintf("𝞧 at time frame=".(k-2))
