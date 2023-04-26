set datafile separator ','
set grid xtics mxtics ytics mytics back
set title 'as^𝛼e^{-bs^𝛽}'
set xlabel 'N'
set ylabel 't'
set key right bottom
f(x) = a*x**b*exp(-c*x**d)

fit f(x) 'PDF_DIAG.csv' via a,b,c,d
set label 1 'model: f(s) = as^𝛼e^{-bs^β}' at 1,0.85 font 'courier,15'
set label 2 sprintf("a=%3.4es", a) at 1,0.80 font 'courier,15'
set label 3 sprintf("𝛼=%3.4es", b) at 1,0.75 font 'courier,15'
set label 4 sprintf("b=%3.4es", c) at 1,0.70 font 'courier,15'
set label 5 sprintf("𝛽=%3.4es", d) at 1,0.65 font 'courier,15'
plot 'PDF_DIAG.csv' w p lc 'red' title 'data', f(x) title 'fitted curve'
