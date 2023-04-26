Read me file Ex5-Lambertini-CODE:

The module with all the subroutine used to complete the different tasks is: Ex_05_module.f90

The main program useful to produce the result is: Ex_05.f90

The two folders "plot_diagonal" and "plot_hermitian" contains all the .csv file produced with the main program and
the gnuplot script useful to build the fits.

R_mid.csv contains the values of <r>.

to compile the module with gfortran compiler:
gfortran -c ex_05_module.f90

to compile the main program with gfortran compiler:
gfortran -o ex05 Ex_05.f90 ex_05_module.o -llapack 
