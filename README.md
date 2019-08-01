# Compile
    $ gcc -std=c99 forcefield?.c bobcode.c  -o bobcode -lfftw3 -lm
# Command-Line Options
    $ ./bobcode
    $ Usage: bobcode dataFile [--nbar relativeCutoff] [--nu accuracyOrder] [-M grid-spacing] 
      [-L numberOfLevels] [--tol-dir direct_tolerance] [--tol-rec reciprocal_tolerance] 
      [--enable-ewald-splitting] [--klim number_of_vawes]

    where data could be any one of *.ini files in the folder without extension.
# Run
    root@www:~/Dropbox/postdoc/src/xcode/bobcode# ./bobcode spceN300
    time_build                     : 0.00238800
    time_energy                    : 0.07388300
    time_total                     : 0.07627100
    data                           : spceN300
    NumberOfLevels                 : 1
    nbar                           : 8.000000
    nu                             : 6
    cutoff                         : 11.428571
    Edge row1                      : 20.00  0.00  0.00
    Edge row2                      :  0.00 20.00  0.00
    Edge row3                      :  0.00  0.00 20.00
    TopLevelMx                     : 7
    TopLevelMy                     : 7
    TopLevelMz                     : 7
    tol_dir                        :  2.441e-05
    tol_rec                        :  1.000e-01
    beta                           : 0.1155843883776609
    kmax                           : 0.0295113638113813
    klimx                          :   2
    klimy                          :   2
    klimz                          :   2
    kLimUserSpecified              :  -1
    effectiveklim_x                :   2
    effectiveklim_y                :   2
    effectiveklim_z                :   2
    ewaldSplitting                 : 0
    utotal                         : -6.4358526791015308e+01
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
    
