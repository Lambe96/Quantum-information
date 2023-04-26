set datafile separator ','
set grid xtics mxtics ytics mytics back
set title 'Fit of time\_loop2.csv'
set xlabel 'N'
set ylabel 't'
set key right bottom
f(x) = a*x**b
a=1e-8
b=3
fit f(x) 'time_loop2.csv' via a,b
set label 1 'model: f(x)= a*x**b' at 200,120 font 'courier,15'
set label 2 sprintf("a=%3.4es", a) at 200,115 font 'courier,15'
set label 3 sprintf("b=%3.4es", b) at 200,110 font 'courier,15'
plot 'time_loop2.csv' w p lc 'red' title 'data', f(x) title 'fitted curve'
