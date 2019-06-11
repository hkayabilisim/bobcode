$ gcc -std=c99 forcefield?.c bobcode.c  -o bobcode -lfftw3 -lm

$ ./bobcode
$ Usage: ./bobcode data nbar nu M L tol_dir tol_rec

$ ./bobcode spceN300 12 6 8 2 0.00001 0.00001
time                           : 1.02310800
data                           : spceN300
L                              : 2
nbar                           : 12.000000
nu                             : 6
cutoff                         : 7.500000
TopLevelMx                     : 8
TopLevelMy                     : 8
TopLevelMz                     : 8
tol_dir                        :  1.000e-05
tol_rec                        :  1.000e-05
beta                           : 0.0916784461234550
kmax                           : 0.0857631311355993
utotal                         : -6.4358736451683683e+01

