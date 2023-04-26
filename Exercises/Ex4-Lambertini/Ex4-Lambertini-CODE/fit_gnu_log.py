#######single-variable fit with gnuplot log-log version###########

import PyGnuplot as gp

data_sep = input('Enter separator character:')
data = input('Enter the name of the file containing the data:')
func = input('Enter the function to fit with gnuplot syntax (f(x)=...):')
num_par = input('Enter the number of parameters:')
num_par = int(num_par)
param_name = []
param_value = []
for i in range(num_par):
    b = input('Enter the name of the parameter #'+str(i+1)+':')
    a = input('Enter the guess for the parameter #'+str(i+1)+':')
    param_value.append(a)
    param_name.append(b)
k=' '
for i in range(num_par):
    k = k+param_name[i]+','
s = list(k)
s = s[:-1]
s = ''.join(s)

tit = input('Enter the title of the plot:')
x_label = input('Enter x label:')
y_label = input('Enter y label:')


gp.c("set datafile separator '"+data_sep+"''")
gp.c('set grid xtics mxtics ytics mytics back ')
gp.c("set title '"+tit+"'")
gp.c('set key right bottom')
gp.c('set xlabel "'+x_label+'"')
gp.c('set ylabel "'+y_label+'"')
gp.c(func)
for i in range(num_par):
    gp.c(param_name[i]+'='+param_value[i])
gp.c("fit f(x) '"+data+"' u (log10($1)):(log10($2)) via"+s)
#gp.c('set label 1 "model: '+func+'" at 2.1,1.5 font "courier,15')
#for i in range(num_par):
    #gp.c('set label '+str(i+2)+'sprintf("'+param_name[i]+'=%3.4f", '+param_name[i]+') at 2.1,'+str(1.3-0.2*i)+' font "courier,15"')
gp.c("plot '"+data+"' u (log10($1)):(log10($2)) w p lc 'red' title 'data', f(x) title 'fitted curve'")
