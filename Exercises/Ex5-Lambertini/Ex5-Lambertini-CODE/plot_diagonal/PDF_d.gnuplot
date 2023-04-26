set datafile separator ','
set grid xtics mxtics ytics mytics back
set title 'Fit of diag\_N\_frac\_10.csv'
set xlabel 's'
set ylabel 'P(s)'
set key right bottom
set xrange [0:10]
f(x) = a*x**b*exp(-c*x**d)

fit f(x) 'diag_N_frac_10.csv' via a,b,c,d
set label 1 'model: f(s) = as^𝛼e^{-bs^β}' at 4,1.2 font 'courier,15'
set label 2 sprintf("a=%3.4es", a) at 4,1.15 font 'courier,15'
set label 3 sprintf("𝛼=%3.4es", b) at 4,1.1 font 'courier,15'
set label 4 sprintf("b=%3.4es", c) at 4,1.05 font 'courier,15'
set label 5 sprintf("𝛽=%3.4es", d) at 4,1 font 'courier,15'
plot 'diag_N_frac_10.csv' w p lc 'red' title 'data', f(x) title 'fitted curve'
