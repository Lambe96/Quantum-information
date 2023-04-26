set datafile separator ','
set grid xtics mxtics ytics mytics back
set title 'Fit of linearized time\_intrinsic.csv'
set xlabel 'log_{10}N'
set ylabel 'log_{10}t'
f(x) = a*x+b
a=3
b=-10
fit f(x) 'time_intrinsic.csv' u (log10($1)):(log10($2)) via a,b
set label 1 'model: f(x)= a*x+b' at 2.2,-0.5 font 'courier,15'
set label 2 sprintf("a=%3.4es", a) at 2.2,-0.7 font 'courier,15'
set label 3 sprintf("b=%3.4es", b) at 2.2,-0.9 font 'courier,15'
plot 'time_intrinsic.csv' u (log10($1)):(log10($2)) w p lc 'red' title 'data', f(x) title 'fitted curve'
