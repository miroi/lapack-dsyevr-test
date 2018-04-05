!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
!
!                            How to compile
!                           ===============
!
! ifort + MKL dynamic lp64:
! ---------------------------
!   ifort -o dsyerv_ifort_mkl.x dsyerv_check.F90  eispack.F -L/mnt/apps/intel/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
!
! gfortran + MKL dynamic lp64 :
! -----------------------------
!   gfortran -o dsyerv_gfortran_mkl.x dsyerv_check.F90  eispack.F -L/mnt/apps/intel/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
!
! pgf90 + MKL dynamic lp64 
! ------------------------
!   pgf90 -o dsyerv_pgf90_mkl.x  dsyerv_check.F90  eispack.F -L/mnt/apps/intel/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
!
! ifort + MKL dynamic ilp64
! -------------------------
!   ifort -o dsyerv_ifort_i8_mkl.x -i8   dsyerv_check.F90  eispack.F -L/mnt/apps/intel/mkl/lib/intel64  -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5
!
! gfortran + MKL dynamic ilp64 : 
! ------------------------------
!  gfortran -o dsyerv_gfortran_i8_mkl.x -fdefault-integer-8  dsyerv_check.F90  eispack.F -L/mnt/apps/intel/mkl/lib/intel64  -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5
!
! ifort + dynamic GNU lp64 libs:
! ------------------------------
!   ifort -o dsyerv_ifort_GNUlibs.x   dsyerv_check.F90  eispack.F  -L/usr/lib64 -lblas -llapack
!
! gfortran + dynamic GNU libs lp64:
! ----------------------------------
!  gfortran -o dsyerv_gfortran_GNUlibs.x  dsyerv_check.F90  eispack.F  -L/usr/lib64 -lblas -llapack
!
! pgf90 +  GNU dynamic lp64 
! ------------------------
!   pgf90 -o dsyerv_pgf90_GNUlibs.x  dsyerv_check.F90 eispack.F  -L/usr/lib64 -lblas -llapack
!
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
module diagmod
!-------------------------------------------------------------------------------------------
!  Module for the standalone testing program of diagonalization routines,
! utils/diag.F
!
! Written by Miro Ilias, August 2014
! Vishnu V. Krishnan, April 2018
!-------------------------------------------------------------------------------------------
  integer :: N=0,NZ=0, print_level
  integer :: lu=12, LUPRI=6
  real*8, allocatable :: Ar(:,:), ASr(:,:),Ai(:,:),ASi(:,:),  & 
                         Aj(:,:),ASj(:,:), Ak(:,:), ASk(:,:), &
                         ArTemp(:,:)
  character*50 :: title_text
  ! hardcoded input matrix name
  character*50 :: matrix_file_name = "Jz_SS_matrix.fermirp2-2"
contains

  subroutine read_matrix
  ! allocates space and reads the matrix for the diagonalization testing from the file
    integer :: i,j, ii,jj
   
    open (lu,file=trim(matrix_file_name), access="sequential", status="old")  
    read(lu,*) N,NZ
    print *,'read N=',N,' NZ=',NZ

    ! allocates space for matrices
    call allocate_space
    ! read matrix elements from its open file
    do i=1,N
    do j=1,N
     if (NZ==1) then
       read(lu,*) ii,jj,Ar(ii,jj)
       ASr(ii,jj)=Ar(ii,jj)
     else if (NZ==2) then
       read(lu,*) ii,jj,Ar(ii,jj),Ai(ii,jj)
       ASr(ii,jj)=Ar(ii,jj); ASi(ii,jj)=Ai(ii,jj)
     else if (NZ==4) then
       read(lu,*) ii,jj,Ar(ii,jj),Ai(ii,jj),Aj(ii,jj),Ak(ii,jj)
       ASr(ii,jj)=Ar(ii,jj); ASi(ii,jj)=Ai(ii,jj)
       ASj(ii,jj)=Aj(ii,jj); ASk(ii,jj)=Ak(ii,jj)
     endif
    enddo
    enddo
    close(lu,status="keep")
    print *,"matrix for diagonalization was read"
 
  end subroutine read_matrix

  subroutine allocate_space
    if (N<=1.or.NZ<=0) then 
       print *,'allocate_space: N=',N,' NZ=',NZ
       stop 5
    endif
    if (NZ==1) then
      allocate( Ar(N,N), stat=ierr )
      allocate( ASr(N,N), stat=ierr )
    else if (NZ==2) then
      allocate( Ar(N,N), Ai(N,N), stat=ierr )
      allocate( ASr(N,N), ASi(N,N), stat=ierr )
    else if (NZ==4) then
      allocate( Ar(N,N), Ai(N,N), Aj(N,N), Ak(N,N), stat=ierr )
      allocate( ASr(N,N), ASi(N,N), ASj(N,N), ASk(N,N), stat=ierr )
    else 
       print *,'wrong value of NZ=',NZ
    endif
  end subroutine allocate_space

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine test_lapack_dsyevr
 !
 ! tests lapack's routine DSYERV for real only symmetric matrices
 !
   implicit none
   integer :: ierr, MATZ, LLWORK, KLWORK, LILWORK
   integer :: i
   integer :: IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, LRA, LRV
   character*1 :: JOBZ, RANGE, UPLO
   real*8 :: ABSTOL, VL, VU
   integer, allocatable :: ISUPPZ(:), IWORK(:)
   real*8, allocatable ::  WORK(:), EIG(:),VEC(:,:)
   real*8 :: DDUM, WORKDUM
   integer :: IWORKDUM,IDUM
   real*8, parameter :: D0=0.0D0

   allocate(ArTemp(N,N), stat=ierr)
   ArTemp = Ar

   ABSTOL=0.0d0; JOBZ='V'; RANGE = 'A'; UPLO = 'U'
   DDUM = D0; IDUM = 0; WORKDUM = D0; IWORKDUM = 0

   MATZ = 1 
   IF(MATZ.EQ.1) THEN
        JOBZ='V'
   ELSEIF(MATZ.EQ.0) THEN
        JOBZ='N'
   ELSE
     WRITE(*,*) 'Illegal value of MATZ: ',MATZ
     STOP 2
   ENDIF

   allocate(EIG(N),VEC(N,N), stat=ierr)
    
   ! 1.call to determine size 
   IERR =  0 
   ! FYI: LLWORK=-1; LIWORK=-1 on input means calculate temporary storage
   M = N ; LRA = N ;LRV = N
   CALL DSYEVR(JOBZ,'A','U',N,ArTemp,LRA,DDUM,DDUM,IDUM,IDUM,ABSTOL,  & 
               M,EIG,VEC,LRV,IDUM,EIG,-1,IWORKDUM,-1,IERR)

   LLWORK = NINT(EIG(1))
   LILWORK = IWORKDUM 

   allocate(WORK(LLWORK), stat=ierr)
   allocate(IWORK(LILWORK), stat=ierr)
   allocate(ISUPPZ(2*N),stat=ierr)

   IERR = 0
   CALL DSYEVR(JOBZ,'A','U',N,ArTemp,LRA,DDUM,DDUM,IDUM,IDUM,ABSTOL,  &
               M,EIG,VEC,LRV,ISUPPZ,WORK,LLWORK,IWORK,LILWORK,IERR)     
   if (IERR /= 0) then
     print *,'The lapack_dsyerv routine ended with error ! ierr=',ierr
   endif

  ! ... print out eigenvalues
  if (print_level >= 0) then
     print *,"LAPACK DSYEVR eigenvalues:"
     do i=1,N
       print *,i,EIG(i)
     enddo
   endif

   call eigv_check(EIG,VEC,'**** LAPACK DSYEVR ****')

   deallocate(ArTemp,WORK,IWORK,EIG,VEC,ISUPPZ,stat=ierr)

  end subroutine test_lapack_dsyevr

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine test_lapack_dstegr
 !
 ! tests lapack's routine DSTEGR for real symmetric matrices
 ! after tri-diagonalising with DSYTRD
 !
   implicit none
   integer :: ierr, MATZ, LLWORK, KLWORK, LILWORK
   integer :: i
   integer :: IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, LRA, LRV
   character*1 :: JOBZ, RANGE, UPLO
   real*8 :: ABSTOL, VL, VU
   integer, allocatable :: ISUPPZ(:), IWORK(:)
   real*8, allocatable ::  WORK(:), EIG(:),VEC(:,:),DIAG(:),OFFDIAG(:), &
       RFLCT(:)
   real*8 :: DDUM, WORKDUM
   integer :: IWORKDUM,IDUM
   real*8, parameter :: D0=0.0D0

   allocate(ArTemp(N,N), stat=ierr)
   ArTemp = Ar

   ABSTOL=0.0d0; JOBZ='V'; RANGE = 'A'; UPLO = 'U'
   DDUM = D0; IDUM = 0; WORKDUM = D0; IWORKDUM = 0

   M = N ; LRA = N ;LRV = N

   allocate(DIAG(N),OFFDIAG(N),RFLCT(N-1))

   IERR =  0 
   CALL DSYTRD(UPLO,N,ArTemp,LRA,DIAG,OFFDIAG(:N-1),RFLCT,WORKDUM,-1,IERR)

   LLWORK = NINT(WORKDUM)
   allocate(WORK(LLWORK))

   IERR = 0
   CALL DSYTRD(UPLO,N,ArTemp,LRA,DIAG,OFFDIAG(:N-1),RFLCT,WORK,LLWORK,IERR)
   if (IERR /= 0) then
     print *,'The lapack_dsytrd routine ended with error ! ierr=',ierr
   endif

   deallocate(RFLCT,WORK)

   MATZ = 1 
   IF(MATZ.EQ.1) THEN
        JOBZ='V'
   ELSEIF(MATZ.EQ.0) THEN
        JOBZ='N'
   ELSE
     WRITE(*,*) 'Illegal value of MATZ: ',MATZ
     STOP 2
   ENDIF

   allocate(EIG(N),VEC(N,N))

   WORKDUM = D0
   ! 1.call to determine size 
   IERR =  0 
   ! FYI: LLWORK=-1; LIWORK=-1 on input means calculate temporary storage
   CALL DSTEGR(JOBZ,'A',N,DIAG,OFFDIAG,DDUM,DDUM,IDUM,IDUM,ABSTOL,  & 
               M,EIG,VEC,LRV,IDUM,WORKDUM,-1,IWORKDUM,-1,IERR)

   LLWORK = NINT(WORKDUM)
   LILWORK = IWORKDUM

   allocate(WORK(LLWORK), stat=ierr)
   allocate(IWORK(LILWORK), stat=ierr)
   allocate(ISUPPZ(2*N),stat=ierr)

   IERR = 0
   CALL DSTEGR(JOBZ,'A',N,DIAG,OFFDIAG,DDUM,DDUM,IDUM,IDUM,ABSTOL,  &
               M,EIG,VEC,LRV,ISUPPZ,WORK,LLWORK,IWORK,LILWORK,IERR)     
   if (IERR /= 0) then
     print *,'The lapack_dstegr routine ended with error ! ierr=',ierr
   endif

  ! ... print out eigenvalues
  if (print_level >= 0) then
     print *,"LAPACK DSTEGR eigenvalues:"
     do i=1,N
       print *,i,EIG(i)
     enddo
   endif

   call eigv_check(EIG,VEC,'**** LAPACK DSYTRD+DSTEGR ****')

   deallocate(ArTemp,DIAG,OFFDIAG,WORK,IWORK,EIG,VEC,ISUPPZ)

  end subroutine test_lapack_dstegr

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  subroutine test_dirac_rs
! tests DIRAC's EISPACK routine RS real Hermitian matrix
    real*8, allocatable :: eigval(:),eigvect(:,:),temp1(:),temp2(:)
    integer :: ierr,i

! ... allocate arrays for the subroutine
    allocate(eigval(N),eigvect(N,N),temp1(N),temp2(N),stat=ierr)
    print *,"arrays for RS subroutine allocated"

    matz=1 ! have both eigenvalues and eigenvectors
    call RS(N,N,Ar,eigval,matz,eigvect,temp1,temp2,ierr)
    if (ierr /= 0) then
      print *,'The RS routine ended with error !'
    endif
! ... print out eigenvalues
    if (print_level >= 0) then
      print *,"Eispack's RS eigenvalues:"
      do i=1,N
        print *,i,eigval(i)
      enddo
    endif

    call eigv_check(eigval,eigvect,'**** EISPACK RS ****')

    deallocate(eigval,eigvect,temp1,temp2,stat=ierr)
  end subroutine test_dirac_rs

  subroutine eigv_check(eigval,eigvect,routine_name)
  !-----------------------------------------------------------------------------------
  !
  ! Do the diagonalization checks using original Hermitian matrix, eigenvalues and
  ! eigenvectors
  ! 
  !           eigvect^+ * AS * eigvect - eigval  = 0
  !
  !           eigvect^+ * eigvect = I
  !
  !           eigvect * eigvect^+ = I
  !
  !-----------------------------------------------------------------------------------
   implicit none
    real*8, intent(in) :: eigval(N),eigvect(N,N)
    character(*), intent(in) :: routine_name
    real*8, allocatable :: AXr(:,:),AYr(:,:),C(:,:)
    integer :: i,j, ierr, iprint
    real*8 :: C_norm_diag, C_norm_offdiag
    real*8 :: AXr_norm_diag, AXr_norm_offdiag
    real*8 :: AYr_norm_diag, AYr_norm_offdiag
    integer(kind=8) :: maxsize = -1
    
    allocate( AXr(N,N),AYr(N,N),C(N,N),stat=ierr )
 
   ! Do  eigvect^{+)*ASr = AXr
    call dgemm('T','N',N,N,N,1.0d0,eigvect,N,ASr,N,0.0d0,AXr,N)

   ! Do  AXr * eigvect = C
    call dgemm('N','N',N,N,N,1.0d0,AXr,N,eigvect,N,0.0d0,C,N)

   ! Another check:  eigvect^{+)*eigvect = I = AXr
    call dgemm('T','N',N,N,N,1.0d0,eigvect,N,eigvect,N,0.0d0,AXr,N)
    if (print_level>=1) then
      write(LUPRI,*) "Printout of eigvect^{+)*eigvect =? I "
    endif

   ! Another check:  eigvect*eigvect^{+} = I = AYr
    call dgemm('N','T',N,N,N,1.0d0,eigvect,N,eigvect,N,0.0d0,AYr,N)
    if (print_level>=1) then
      write(LUPRI,*) "Printout of eigvect*eigvect^{+} =? I "
    endif

    ! substract eigenvalues ...
    C_norm_diag = 0.0d0; AXr_norm_diag = 0.0d0; AYr_norm_diag = 0.0d0
    do i=1,N
      C(i,i) = C(i,i) -eigval(i); 
      AXr(i,i)=AXr(i,i)-1.0d0; AYr(i,i)=AYr(i,i)-1.0d0
      C_norm_diag = C_norm_diag + DABS(C(i,i))
      AXr_norm_diag = AXr_norm_diag  + DABS(AXr(i,i))
      AYr_norm_diag = AYr_norm_diag  + DABS(AYr(i,i))
    enddo
    C_norm_diag   = C_norm_diag/DFLOAT(N)
    AXr_norm_diag = AXr_norm_diag/DFLOAT(N)
    AYr_norm_diag = AYr_norm_diag/DFLOAT(N)

    C_norm_offdiag = 0.0d0; AXr_norm_offdiag = 0.0d0;  AYr_norm_offdiag = 0.0d0
    do i = 1, N
    do j = 1, N
     if (i /= j) then
       C_norm_offdiag = C_norm_offdiag + DABS(C(i,j))
       AXr_norm_offdiag = AXr_norm_offdiag + DABS(AXr(i,j))
       AYr_norm_offdiag = AYr_norm_offdiag + DABS(AYr(i,j))
     endif
    enddo
    enddo
    C_norm_offdiag = C_norm_offdiag/((DFLOAT(N)*DFLOAT(N))-DFLOAT(N))
    AXr_norm_offdiag = AXr_norm_offdiag/((DFLOAT(N)*DFLOAT(N))-DFLOAT(N))
    AYr_norm_offdiag = AYr_norm_offdiag/((DFLOAT(N)*DFLOAT(N))-DFLOAT(N))

    write(LUPRI,"(2X,A)") routine_name
    write(LUPRI,"(1X,A,D10.4,2X,A,D10.4)") & 
 "U^{+}*A*U - eps ?= 0> norm/diag:",C_norm_diag,  "norm/offdiag:",C_norm_offdiag

    write(LUPRI,"(1X,A,D10.4,2X,A,D10.4)") & 
 "    U^{+}*U - I ?= 0> norm/diag:",AXr_norm_diag,"norm/offdiag:",AXr_norm_offdiag
    write(LUPRI,"(1X,A,D10.4,2X,A,D10.4)") & 
 "    U*U^{+} - I ?= 0> norm/diag:",AYr_norm_diag,"norm/offdiag:",  &
  AYr_norm_offdiag


    deallocate( AXr, AYr, C, stat=ierr )

  end subroutine eigv_check

end module diagmod


Program Test_DIRAC_Diagonalization_Routines
  use diagmod
   write(6,'(2x,a)') "Welcome to the testing program for DIRAC's diagonalization routines !"

  call read_matrix
  call test_dirac_rs
  call test_lapack_dsyevr
  call test_lapack_dstegr

End Program Test_DIRAC_Diagonalization_Routines
