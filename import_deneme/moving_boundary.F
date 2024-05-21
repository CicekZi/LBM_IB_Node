        SUBROUTINE MOVING_BC()

#ifdef MOVING_BOUNDARY


        IMPLICIT NONE

        ! Define precision
        INTEGER, PARAMETER :: WP = KIND(1.0D0)

        ! Allocatable arrays
        REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: RBLANK
        REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: ZC, D, ZC_MIN, DELTA
        REAL(WP), ALLOCATABLE, DIMENSION(:)     :: XX, ZZ
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IJK
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: NUP_TOT
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NUP

        ! Scalar variables
        REAL(WP) :: RAD, X0, X1, X2, XW, ZW, XB, ZB, XINT
        REAL(WP) :: X00, X11, X22
        REAL(WP) :: CENTER, TANDEM_CYL, DIFFER
        INTEGER  :: I, J, K, IGL, JGL, KGL, NC, NCURV, NCURV_BUF, NADV_DEN, TOP
        INTEGER :: io_status

        ! File names
        CHARACTER(LEN=20) :: file_NUP, file_NUP_TOT, file_D, file_DELTA 
        CHARACTER(LEN=20) :: file_RBLANK, file_IJK, file_NCURV

        ! Assign file names
        file_NUP     = 'NUP.txt'
        file_NUP_TOT = 'NUP_TOT.txt'
        file_D       = 'D.txt'
        file_DELTA   = 'DELTA.txt'
        file_RBLANK  = 'RBLANK.txt'
        file_IJK     = 'IJK.txt'
        file_NCURV   = 'NCURV.txt'

   
        

        ! Allocate arrays with fixed sizes
        ALLOCATE (IJK(1:3,1:NCURV_BUF))
        ALLOCATE (NUP_TOT(1:NCURV_BUF))
        ALLOCATE (NUP(1:18,1:NCURV_BUF))
        ALLOCATE (D(1:18,1:NCURV))
        ALLOCATE (DELTA(1:18,1:NCURV))
        ALLOCATE (RBLANK(-2:303, -2:303, -2:203))

        ! Initialize arrays to zero
        IJK     = 0
        NUP_TOT = 0
        NUP     = 0
        D = 0.0_WP
        DELTA   = 0.0_WP
        RBLANK = 0.0_WP

        ! Read data from files
        open(unit=10, file=file_NCURV, status='old')
        if (io_status /= 0) then
            print *, "Error opening ", file_NCURV
            stop
        end if
        read(10, *) NCURV
        close(10)

        open(unit=20, file=file_NUP, status='old')
        if (io_status /= 0) then
            print *, "Error opening ", file_NUP
            stop
        end if
        read(20, *) NUP
        close(20)

        open(unit=30, file=file_NUP_TOT, status='old')
        if (io_status /= 0) then
            print *, "Error opening ", file_NUP_TOT
            stop
        end if
        read(30, *) NUP_TOT
        close(30)

        open(unit=40, file=file_D, status='old')
        if (io_status /= 0) then
            print *, "Error opening ", file_D
            stop
        end if
        read(40, *) D
        close(40)

        open(unit=50, file=file_DELTA, status='old')
        if (io_status /= 0) then
            print *, "Error opening ", file_DELTA
            stop
        end if
        read(50, *) DELTA
        close(50)

        open(unit=60, file=file_RBLANK, status='old')
        if (io_status /= 0) then
            print *, "Error opening ", file_RBLANK
            stop
        end if
        read(60, *) RBLANK
        close(60)

        open(unit=70, file=file_IJK, status='old')
        if (io_status /= 0) then
            print *, "Error opening ", file_IJK
            stop
        end if
        read(70, *) IJK
        close(70)

#endif

        RETURN
        END SUBROUTINE MOVING_BC
