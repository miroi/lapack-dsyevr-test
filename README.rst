==================
lapack-dsyevr-test
==================

Important self-test of the famous lapack's diagonalization routine DSYEVR. It checks eigenvectors orthogonality. GNU lapack is affected with the bug, see https://github.com/Reference-LAPACK/lapack/issues/151. 

This small software buildup is based upon autocmake project ( https://github.com/coderefinery/autocmake ).

Quick buildup
-------------

Clone the repository, go the src/ directory and type there:

::

 gfortran -fallow-argument-mismatch dsyevr_check.F90 eispack.F -llapack -lblas

you get the a.out executable you can work with.

Buildup examples 
----------------

- Gfortran 4.7.2 and later versions, with system native libraries (/usr/lib/liblapack.so /usr/lib/libblas.so)

::

 python setup --fc=gfortran --extra-fc-flags="-fallow-argument-mismatch" --cmake-options="-D MATH_LIB_SEARCH_ORDER='SYSTEM_NATIVE'"  build_gfortran_sysnatlibs
 cd build_gfortran_sysnatlibs
 make VERBOSE=1
 cd src
 ./dsyerv_check

confirms eigenvectors nonorthogonality in the dsyevr routine:

::

 **** LAPACK DSYEVR ****
 U^{+}*A*U - eps ?= 0> norm/diag:0.1147D-14  norm/offdiag:0.1365D-05
     U^{+}*U - I ?= 0> norm/diag:0.3824D-15  norm/offdiag:0.2729D-05
     U*U^{+} - I ?= 0> norm/diag:0.8556D-05  norm/offdiag:0.5285D-05


- Gfortran (4.7.2) and the own downloaded Lapack 3.6.0

::
 
 python setup --fc=gfortran --blas=off --lapack=off --cmake-options="-D EXPLICIT_LIBS='-L/u/milias/Work/qch/software/lapack/lapack-3.6.0/build/lib -llapack -lblas'"  build_gfortran_lapack3.6.0
 cd build_gfortran_lapack3.6.0
 make VERBOSE=1
 cd src
 ./dsyerv_check


shows the nonorthogonality:

::

  **** LAPACK DSYEVR ****
 U^{+}*A*U - eps ?= 0> norm/diag:0.1616D-14  norm/offdiag:0.1419D-06
     U^{+}*U - I ?= 0> norm/diag:0.3454D-15  norm/offdiag:0.2838D-06
     U*U^{+} - I ?= 0> norm/diag:0.1047D-05  norm/offdiag:0.7950D-06


- Gfortran (4.4.7) with the fresh Lapack 3.7.0

::

 python setup --fc=gfortran --blas=off --lapack=off --cmake-options="-D EXPLICIT_LIBS='-L/home/milias/Work/qch/software/smaller_software_projects/lapack-dsyevr-test/lapack-3.7.0/build/lib  -llapack -lblas'"  build_gfortran_lapack3.7.0
 cd build_gfortran_lapack3.7.0
 make VERBOSE=1
 cd src
 ./dsyerv_check

shows the nonorthogonality:

::

  **** LAPACK DSYEVR ****
 U^{+}*A*U - eps ?= 0> norm/diag:0.1616D-14  norm/offdiag:0.1419D-06
     U^{+}*U - I ?= 0> norm/diag:0.3454D-15  norm/offdiag:0.2838D-06
     U*U^{+} - I ?= 0> norm/diag:0.1047D-05  norm/offdiag:0.7950D-06


- Intel 15.0 and the downloaded Lapack 3.6.0

::

 python setup --fc=ifort --blas=off --lapack=off --cmake-options="-D EXPLICIT_LIBS='-L/u/milias/Work/qch/software/lapack/lapack-3.6.0/build/lib -llapack -lblas -lgfortran'"  build_ifort_lapack3.6.0
 cd build_ifort_lapack3.6.0
 make VERBOSE=1
 cd src
 ./dsyerv_check

confirms eigenvectors orthogonality in the dsyevr routine:

::

  **** LAPACK DSYEVR ****
 U^{+}*A*U - eps ?= 0> norm/diag:0.1604D-14  norm/offdiag:0.2574D-15
     U^{+}*U - I ?= 0> norm/diag:0.2961D-15  norm/offdiag:0.9328D-16
     U*U^{+} - I ?= 0> norm/diag:0.1974D-15  norm/offdiag:0.1216D-15


- Intel 15.0/16.0 with the attached MKL library

::

 python setup --fc=ifort --cmake-options="-D MATH_LIB_SEARCH_ORDER='MKL'"  build_ifort_mkl
 cd build_ifort_mkl
 make VERBOSE=1
 cd src
 ./dsyerv_check

has proper orthogonality in the dsyevr routine as we can see from the output:

::
 
   **** LAPACK DSYEVR ****
 U^{+}*A*U - eps ?= 0> norm/diag:0.1684D-14  norm/offdiag:0.2746D-15
     U^{+}*U - I ?= 0> norm/diag:0.2097D-15  norm/offdiag:0.6923D-16
     U*U^{+} - I ?= 0> norm/diag:0.2467D-16  norm/offdiag:0.8681D-16

