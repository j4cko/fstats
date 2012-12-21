fstats
======

calculates statistical quantities on the commandline, currently including mean, naive sd, bootstrapped sd, binning analysis (to be continued). can be used stand-alone or as Fortran module

FEATURES
========

- arithmetic mean
- empirical standard deviation / error
- autocorrelation function
- integrated autocorrelation time
- binning: error estimation via emp. se or bootstrap

BUILD
=====

use your favorite fortran compiler, e.g. gfortran (>4.6 or something like that)

gfortran -c fstats_mod.f90
gfortran -c fstats.f90
gfortran fstats_mod.o fstats.o -o fstats

USE
===

at the moment I recommend use as a module, look in fstats_mod.f90 to see which routines you can use.

TODO
====

- make frontend more flexible, commandline options
- maybe: text/graphical output of autocorr. fun
- documentation
