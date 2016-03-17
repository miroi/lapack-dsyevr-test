==================
lapack-dsyevr-test
==================

Important self-test of the famous lapack's diagonalization routine DSYEVR.

Software buildup is based upon autocmake project ( https://github.com/scisoft/autocmake ).

Buildup examples 
----------------

- Gfortran and system native libraries (/usr/lib/liblapack.so /usr/lib/libblas.so)

::

 python setup --fc=gfortran --cmake-options="-D MATH_LIB_SEARCH_ORDER='SYSTEM_NATIVE'"  build_gfortran_sysnatlibs
 cd build_gfortran_sysnatlibs
 make VERBOSE=1
 cd src
 ./dsyerv_check

confirms uncorrect dsyevr routine:

::

 **** LAPACK DSYEVR ****
 U^{+}*A*U - eps ?= 0> norm/diag:0.1147D-14  norm/offdiag:0.1365D-05
     U^{+}*U - I ?= 0> norm/diag:0.3824D-15  norm/offdiag:0.2729D-05
     U*U^{+} - I ?= 0> norm/diag:0.8556D-05  norm/offdiag:0.5285D-05

- Gfortran and downloaded Lapack 3.6.0

::
 
 python setup --fc=gfortran --blas=off --lapack=off --cmake-options="-D EXPLICIT_LIBS='-L/u/milias/Work/qch/software/lapack/lapack-3.6.0/build/lib -llapack -lblas'"  build_gfortran_lapack3.6.0
 cd build_gfortran_lapack3.6.0
 make VERBOSE=1
 cd src
 ./dsyerv_check


shows wrong dsyevr routine:

::

  **** LAPACK DSYEVR ****
 U^{+}*A*U - eps ?= 0> norm/diag:0.1616D-14  norm/offdiag:0.1419D-06
     U^{+}*U - I ?= 0> norm/diag:0.3454D-15  norm/offdiag:0.2838D-06
     U*U^{+} - I ?= 0> norm/diag:0.1047D-05  norm/offdiag:0.7950D-06

- Intel 15.0/16.0 with MKL

::

 python setup --fc=ifort --cmake-options="-D MATH_LIB_SEARCH_ORDER='MKL'"  build_ifort_mkl
 cd build_ifort_mkl
 make VERBOSE=1
 cd src
 ./dsyerv_check

has proper dsyevr routine as we can see from the output

::
 
   **** LAPACK DSYEVR ****
 U^{+}*A*U - eps ?= 0> norm/diag:0.1684D-14  norm/offdiag:0.2746D-15
     U^{+}*U - I ?= 0> norm/diag:0.2097D-15  norm/offdiag:0.6923D-16
     U*U^{+} - I ?= 0> norm/diag:0.2467D-16  norm/offdiag:0.8681D-16

