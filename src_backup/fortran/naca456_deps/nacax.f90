!+
MODULE NacaAuxiliary
! ------------------------------------------------------------------------------
! PURPOSE -
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!
! REVISION HISTORY
!   DATE  VERS RERSON  STATEMENT OF CHANGES
! 17Jun99  0.5   RLC   Original coding
!  6Sep99  0.6   RLC   Added LeadingEdgeRadius function
!  9Dec99  0.7   RLC   Moved all epsilon, psi tables to new module
! 13Feb01  0.8   RLC   Added CombineThicknessAndCamber
! 16Feb01  0.81  RLC   Added ParametrizeAirfoil, Combo6seriesMeanLine
! 14Mar01  0.85  RLC   Added Polynomial,ScaleFactor;revised SetSixDigitPoints
! 15Mar01  0.86  RLC   Added spline interpolation to 6-series thickness
! 19Sep01  0.9   RLC   Revised LoadX - coarse is same as Abbot & vonDoenhoeff
! 21Sep01  0.91  RLC   Added MEDIUM and FINE spacing to LoadX
! 25Sep01  0.92  RLC   Added GetRk1 and GetRk1k2
! 03Oct01  0.93  RLC   Corrected two errors in 3-digit reflex mean line
! 13Oct01  0.94  RLC   Revised CalA1,CalA2,CalA3 
! 16Nov01  0.95  RLC   Getting the interpolation right
! 24Nov01  0.96  RLC   Optional argument coding
! 30Nov01  0.97  RLC   Added LeRadius4 and LeRadius4M
! 06Dec01  0.98  RLC   Added LeRadius6
! 26Dec01  0.99  RLC   Made sure we could do points dense at l.e.
!  4Jan02  1.00  RLC   Final cleanup for PDAS 7
! 16Jan09  1.05  RLC   Small cosmetic improvements
! 07Nov10  1.10  RLC   Fixed bug in GetRk1 (thanks to Robert Stone)

IMPLICIT NONE

  CHARACTER(LEN=*),PUBLIC,PARAMETER:: AUX_VERSION = "1.10 (7 November 2010)"

  PRIVATE:: AddTrailingEdgePointIfNeeded
  PRIVATE:: CalA1,CalA2,CalA3
  PRIVATE:: CalculateD1
  PUBLIC:: CombineThicknessAndCamber
  PUBLIC:: GetRk1             ! used by avd
  PUBLIC:: GetRk1K2           ! used by avd
  PUBLIC:: InterpolateCombinedAirfoil
  PUBLIC:: InterpolateUpperAndLower
  PUBLIC:: LeadingEdgeRadius4
  PUBLIC:: LeadingEdgeRadius4M
  PUBLIC:: LeadingEdgeRadius6
  PUBLIC:: MeanLine2
  PUBLIC:: MeanLine3
  PUBLIC:: MeanLine3reflex
  PUBLIC:: MeanLine6
  PUBLIC:: MeanLine6M
  PUBLIC:: ParametrizeAirfoil
  PRIVATE:: Polynomial
  PUBLIC:: ScaleFactor
  PUBLIC:: SetSixDigitPoints
  PUBLIC:: Thickness4
  PUBLIC:: Thickness4M
  PUBLIC:: Thickness4sharpTE
  PUBLIC:: Thickness6


CONTAINS

!+
SUBROUTINE AddTrailingEdgePointIfNeeded(n,x,y)
! ------------------------------------------------------------------------------
! PURPOSE - If x(1:n) and y(1:n) define an airfoil surface from leading
!   edge to trailing edge and if the final point has x < 1.0, add an
!   additional point that extrapolates to x=1.0
!   x and y must be dimensioned at least n+1 or no action is taken
  INTEGER,INTENT(IN OUT):: n
  REAL,INTENT(IN OUT),DIMENSION(:):: x,y

  REAL:: slope
!-------------------------------------------------------------------------------
  IF (SIZE(x) <= n) RETURN
  IF (SIZE(y) <= n) RETURN
  IF (x(n) >= 1.0) RETURN

  slope=(y(n)-y(n-1))/(x(n)-x(n-1))
  x(n+1)=1.0
  y(n+1)=y(n)+slope*(1.0-x(n))
  n=n+1
  RETURN
END Subroutine AddTrailingEdgePointIfNeeded   ! --------------------------------

!+
PURE FUNCTION CalA1(a0, a2, a3, xmt) RESULT(f)
! ------------------------------------------------------------------------------
!  PURPOSE - The a1 variable used for the modified 4-digit thickness   
  REAL,INTENT(IN):: a0,a2,a3,xmt
  REAL:: f
  REAL:: v1,v2,v3
!-------------------------------------------------------------------------------
  v1=0.5*a0/SQRT(xmt)
  v2=2.0*a2*xmt 
  v3=3.0*a3*xmt*xmt
  f= -v1 - v2 - v3 
  RETURN 
END Function CalA1   ! ---------------------------------------------------------
                                                                        
!+
PURE FUNCTION CalA2(a0, a3, xmt) RESULT(f)
! ------------------------------------------------------------------------------
! PURPOSE - The a2 variable used for the modified 4-digit thickness   
  REAL,INTENT(IN):: a0,a3,xmt

  REAL:: f
  REAL:: v1,v2,v3
!-------------------------------------------------------------------------------
  v1=0.1/(xmt*xmt)
  v2=0.5*a0/SQRT(xmt*xmt*xmt)
  v3=2.0*a3*xmt                                                                        
  f = -v1 + v2 - v3 
  RETURN 
END Function Cala2   ! ---------------------------------------------------------

!+
PURE FUNCTION CalA3(a0, d1, xmt) RESULT(f)
! ------------------------------------------------------------------------------
! PURPOSE - The a3 variable used for the modified 4-digit thickness   

  REAL,INTENT(IN):: a0,d1,xmt

  REAL:: f
  REAL:: omxmt
  REAL:: v1,v2,v3                                                   
!-------------------------------------------------------------------------------
  omxmt= 1.0-xmt 
  v1=0.1/(xmt*xmt*xmt)
  v2=(d1*omxmt-0.294)/(xmt*omxmt*omxmt) 
  v3=(3.0/8.0)*a0/(xmt**2.5) 
  f = v1 + v2 - v3 
  RETURN
END Function CalA3   ! ---------------------------------------------------------

!+
PURE FUNCTION CalculateD1(m) RESULT(f)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate d1, the trailing edge half angle, used for the
!   modified 4-digit thickness distribution.
! METHOD - Curve fit of fourth order polynomial to the data:
!   m     d1
!  0.2  0.200
!  0.3  0.234
!  0.4  0.315
!  0.5  0.465
!  0.6  0.700

  REAL,INTENT(IN):: m  ! x-location of max thickness, fraction of chord
  REAL:: f  ! d1

  REAL,PARAMETER,DIMENSION(5):: A= (/ &   ! A(1)=constant;A(2)=linear,etc.
    3.48E-5, 2.3076628, -10.127712, 19.961478, -10.420597 /)
!-------------------------------------------------------------------------------
! OLD METHOD: (not too good)
!  REAL,PARAMETER:: C1= 2.24, C2= 5.42, C3= 12.3, C4= 0.878
!  xmtsq= xmt * xmt
!  f = 0.1 * (c1 - (c2*xmt) + (c3*xmtsq)) / (1.0 - c4*xmt)
 
  f=Polynomial(A,m)

  RETURN
END Function CalculateD1   ! ---------------------------------------------------

!+
SUBROUTINE CombineThicknessAndCamber(x,thick,y,yp, &
  xupper,yupper, xlower,ylower)
! NOTE - The dimension of xupper,... must not be less than that of x.
! ------------------------------------------------------------------------------
! PURPOSE - Add the computed thickness and camber to get a cambered airfoil.

  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(IN),DIMENSION(:):: thick   ! note this the half-thickness
  REAL,INTENT(IN),DIMENSION(:):: y
  REAL,INTENT(IN),DIMENSION(:):: yp

  REAL,INTENT(OUT),DIMENSION(:):: xupper
  REAL,INTENT(OUT),DIMENSION(:):: yupper
  REAL,INTENT(OUT),DIMENSION(:):: xlower
  REAL,INTENT(OUT),DIMENSION(:):: ylower

  INTEGER:: n
  REAL,ALLOCATABLE,DIMENSION(:):: s,c
!-------------------------------------------------------------------------------
  n=SIZE(x)
  ALLOCATE(s(n), c(n))
  s=SIN(ATAN(yp))
  c=COS(ATAN(yp))

  xupper(1:n)=x-thick*s   ! will fail if dimension of xupper is too small
  yupper(1:n)=y+thick*c

  xlower(1:n)=x+thick*s
  ylower(1:n)=y-thick*c

  DEALLOCATE(c,s)
  RETURN
END Subroutine CombineThicknessAndCamber   ! -----------------------------------

!+
SUBROUTINE GetRk1(x,r,k1)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the r and k1 factors used in the calculation of
!   ordinates of a three-digit mean line
! REF: Table, upper left of p.8 of NASA Technical Memorandum 4741
USE SplineProcedures,ONLY: TableLookup  
  REAL,INTENT(IN):: x   ! x-coor of max. camber ( in [0,0.25] )
  REAL,INTENT(OUT):: r   ! factor used by MeanLine3
  REAL,INTENT(OUT):: k1  ! factor used by MeanLine3

  REAL,PARAMETER,DIMENSION(5)::    M = (/  0.05,   0.1,   0.15,   0.2, 0.25/)
  REAL,PARAMETER,DIMENSION(5):: RTAB = (/0.0580, 0.126, 0.2025,  0.29,0.391/)
  REAL,PARAMETER,DIMENSION(5):: KTAB = (/ 361.4, 51.64, 15.957, 6.643, 3.23/)
!                  typo fixed 7 Nov 2010   -----------------^ 
!                  thanks to Robert Stone 
!-------------------------------------------------------------------------------
  r=TableLookup(M,RTAB,1,x)
  k1=TableLookup(M,KTAB,1,x)
  RETURN
END Subroutine GetRk1   ! ------------------------------------------------------

!+
SUBROUTINE GetRk1k2(x,r,k1,k2)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the r, k1 and k2 factors used in the calculation of
!   ordinates of a three-digit-reflex mean line
! REF - Table, upper right of p.8 of NASA Technical Memorandum 4741.
! NOTE - This is really the correct equation for k2. It was printed
!   incorrectly on p. 9 of NASA TM X-3284 in 1976 and then 
!   changed (incorrectly!!) on p.8 of NASA TM 4741 in 1995.  Eastman Jacobs 
!   had it right the first time on p.522 of NACA Report 537 in 1935.
!*** USE SplineProcedures,ONLY: TableLookup
USE SplineProcedures,ONLY: TableLookup
  REAL,INTENT(IN):: x   ! x-coor of maximum camber  ( in [0,0.25] )
  REAL,INTENT(OUT):: r  ! factor used by MeanLine3R
  REAL,INTENT(OUT):: k1 ! factor used by MeanLine3R
  REAL,INTENT(OUT):: k2 ! factor used by MeanLine3R

  REAL,PARAMETER,DIMENSION(4):: M = (/ 0.1,0.15,0.2,0.25 /)
  REAL,PARAMETER,DIMENSION(4):: RTAB = (/0.13,0.217,0.318,0.441 /)
  REAL,PARAMETER,DIMENSION(4):: KTAB = (/ 51.99,15.793,6.52,3.191 /)
!  REAL,PARAMETER,DIMENSION(4):: K2TAB = (/ 0.000764,0.00677,0.0303,0.1355 /)
!-------------------------------------------------------------------------------
  r=TableLookup(M,RTAB,1,x)
  k1=TableLookup(M,KTAB,1,x)
  k2=(3.0*(r-x)**2-r**3)/(1.0-r)**3
  RETURN
END Subroutine GetRk1k2   ! ----------------------------------------------------

!+
SUBROUTINE InterpolateCombinedAirfoil(x,yt,ymean,ymeanp, yu,yl, yup,ylp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the upper and lower cordinates of an airfoil defined by
!   the thickness (x,yt) and the mean line (x,ymean,ymeanp). 
USE SplineProcedures,ONLY: FMMspline,SplineZero,PClookup
  REAL,INTENT(IN),DIMENSION(:):: x,yt,ymean,ymeanp
  REAL,INTENT(OUT),DIMENSION(:):: yu,yl
  REAL,INTENT(OUT),DIMENSION(:),OPTIONAL:: yup,ylp

  REAL:: dxds,dyds
!  REAL:: dummy
  INTEGER:: errCode
  INTEGER:: k
!  REAL:: maxError
  INTEGER:: n,nn
  INTEGER:: nupper,nlower
  REAL:: sbar
  REAL,PARAMETER:: TOL=1E-6
  REAL,ALLOCATABLE,DIMENSION(:):: xupper,yupper,xlower,ylower
  REAL,ALLOCATABLE,DIMENSION(:):: xLocal,yLocal,sLocal,xpLocal,ypLocal
!-------------------------------------------------------------------------------
  n=SIZE(x)   

  ALLOCATE(xupper(n+1),yupper(n+1),xlower(n+1),ylower(n+1))
  CALL CombineThicknessAndCamber(x,yt,ymean,ymeanp, &
     xupper,yupper,xlower,ylower)
  nupper=n
  CALL AddTrailingEdgePointIfNeeded(nupper,xupper,yupper) ! might add 1 point
  nlower=n
  CALL AddTrailingEdgePointIfNeeded(nlower,xlower,ylower) ! might add 1 point

  nn=nupper+nlower-1
  ALLOCATE(xLocal(nn),yLocal(nn),sLocal(nn),xpLocal(nn),ypLocal(nn))
  CALL ParametrizeAirfoil(xupper(1:nupper),yupper(1:nupper), &
    xlower(1:nlower),ylower(1:nlower), sLocal,xLocal,yLocal)
  DEALLOCATE(xupper,yupper,xlower,ylower)
!!!  WRITE(*,*) "Airfoil parametrized"   ! omit after checkout
  
!... now fit splines to xLocal vs. sLocal  and  yLocal vs. slocal
!.   This is a total fit around upper and lower surfaces
  CALL FMMSpline(sLocal,xLocal,xpLocal)
  CALL FMMSpline(sLocal,yLocal,ypLocal)
!  WRITE(DBG,*) "total airfoil parameterized data"
!  WRITE(DBG,*) "         sLocal    xLocal    xpLocal    yLocal    ypLocal"
!  CALL PrintArraysNumbered(DBG, sLocal, xLocal,xpLocal, yLocal,ypLocal)
  
!  nx=SIZE(x)
!!!!  ALLOCATE(yu(nx),yl(nx),yup(nx),ylp(nx))

!  maxError=0.0
  DO k=1,n
    CALL SplineZero(sLocal(1:nupper),xLocal(1:nupper),xpLocal(1:nupper), &
      x(k),TOL,sbar,errCode)
    IF (errCode /= 0) THEN
      WRITE(*,*) "errCode not zero from SplineZero"
    END IF
    CALL PClookup(slocal(1:nupper),xLocal(1:nupper),xpLocal(1:nupper), &
      sbar,FP=dxds)
    CALL PClookup(slocal(1:nupper),yLocal(1:nupper),ypLocal(1:nupper), &
      sbar,yu(k),dyds)
    IF (Present(yup)) THEN
      IF (dxds==0.0) THEN
        yup(k)=0.0
      ELSE
        yup(k)=dyds/dxds
      END IF
    END IF

    CALL SplineZero(sLocal(nupper:nn),xLocal(nupper:nn),xpLocal(nupper:nn), &
      x(k),TOL,sbar,errCode)
    IF (errCode /= 0) THEN
      WRITE(*,*) "errCode not zero"
    END IF
    CALL PClookup(slocal(nupper:nn),xLocal(nupper:nn),xpLocal(nupper:nn), &
      sbar,FP=dxds)
    CALL PClookup(slocal(nupper:nn),yLocal(nupper:nn),ypLocal(nupper:nn), &
      sbar,yl(k),dyds)
    IF (Present(ylp)) THEN
      IF (dxds==0.0) THEN
        ylp(k)=0.0
      ELSE
        ylp(k)=dyds/dxds
      END IF
    END IF
  END DO

  DEALLOCATE(ypLocal,xpLocal,sLocal,yLocal,xLocal)

  RETURN
END Subroutine InterpolateCombinedAirfoil   ! ----------------------------------

!+
SUBROUTINE InterpolateUpperAndLower(xupper,yupper, xlower,ylower,       & 
  x,yu,yl, yup,ylp)
! ------------------------------------------------------------------------------
! PURPOSE - Using the airfoil defined by (xupper,yupper) for upper surface
!   and (xlower,ylower) for the lower surface, interpolate at each point of
!   array x to yield yu and yup on upper surface and yl and ylp on lower.
! NOTE - xupper and xlower do not need to be the same size.
USE SplineProcedures,ONLY: FMMspline,SplineZero,PClookup
  REAL,INTENT(IN),DIMENSION(:):: xupper,yupper
  REAL,INTENT(IN),DIMENSION(:):: xlower,ylower
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(OUT),DIMENSION(:):: yu,yl
  REAL,INTENT(OUT),DIMENSION(:),OPTIONAL:: yup,ylp

  REAL:: dxds,dyds
!  REAL:: dummy
  INTEGER:: errCode
  INTEGER:: k
!  REAL:: maxError
  INTEGER:: nn
  INTEGER:: nupper,nlower
  REAL:: sbar
  REAL,PARAMETER:: TOL=1E-6

  REAL,ALLOCATABLE,DIMENSION(:):: xupperCopy,yupperCopy,xlowerCopy,ylowerCopy
  REAL,ALLOCATABLE,DIMENSION(:):: xLocal,yLocal,sLocal,xpLocal,ypLocal
!-------------------------------------------------------------------------------
  nupper=SIZE(xupper)
  nlower=SIZE(xlower)
  ALLOCATE(xupperCopy(nupper+1), yupperCopy(nupper+1))
  ALLOCATE(xlowerCopy(nlower+1), ylowerCopy(nlower+1))
  xupperCopy(1:nupper)=xupper(1:nupper)
  yupperCopy(1:nupper)=yupper(1:nupper)
  xlowerCopy(1:nlower)=xlower(1:nlower)
  ylowerCopy(1:nlower)=ylower(1:nlower)
  CALL AddTrailingEdgePointIfNeeded(nupper,xupperCopy,yupperCopy)
  CALL AddTrailingEdgePointIfNeeded(nlower,xlowerCopy,ylowerCopy)

  nn=nupper+nlower-1
  
  ALLOCATE(xLocal(nn),yLocal(nn),sLocal(nn),xpLocal(nn),ypLocal(nn))
  CALL ParametrizeAirfoil(xupperCopy(1:nupper),yupperCopy(1:nupper), &
    xlowerCopy(1:nlower),ylowerCopy(1:nlower), sLocal,xLocal,yLocal)
  DEALLOCATE(ylowerCopy,xlowerCopy,yupperCopy,xupperCopy)
  
!... now fit splines to xLocal vs. sLocal  and  yLocal vs. slocal
!.   This is a total fit around upper and lower surfaces
  CALL FMMSpline(sLocal,xLocal,xpLocal)
  CALL FMMSpline(sLocal,yLocal,ypLocal)

!  maxError=0.0
  DO k=1,SIZE(x)
    CALL SplineZero(sLocal(1:nupper),xLocal(1:nupper),xpLocal(1:nupper), &
      x(k),TOL,sbar,errCode)
    IF (errCode /= 0) THEN
      WRITE(*,*) "errCode not zero"
    END IF
    CALL PClookup(slocal(1:nupper),xLocal(1:nupper),xpLocal(1:nupper), &
      sbar,FP=dxds)
    CALL PClookup(slocal(1:nupper),yLocal(1:nupper),ypLocal(1:nupper), &
      sbar,yu(k),dyds)
    IF (Present(yup)) THEN
      IF (dxds==0.0) THEN
        yup(k)=0.0
      ELSE
        yup(k)=dyds/dxds
      END IF
    END IF

    CALL SplineZero(sLocal(nupper:nn),xLocal(nupper:nn),xpLocal(nupper:nn), &
      x(k),TOL,sbar,errCode)
    IF (errCode /= 0) THEN
      WRITE(*,*) "errCode not zero"
    END IF
    CALL PClookup(slocal(nupper:nn),xLocal(nupper:nn),xpLocal(nupper:nn), &
      sbar,FP=dxds)
    CALL PClookup(slocal(nupper:nn),yLocal(nupper:nn),ypLocal(nupper:nn), &
      sbar,yl(k),dyds)
    IF (Present(ylp)) THEN
      IF (dxds==0.0) THEN
        ylp(k)=0.0
      ELSE
        ylp(k)=dyds/dxds
      END IF
    END IF
  END DO

  DEALLOCATE(ypLocal,xpLocal,sLocal,yLocal,xLocal)

  RETURN
END Subroutine InterpolateUpperAndLower   ! ------------------------------------

!+
PURE FUNCTION LeadingEdgeRadius4(toc) RESULT(rle)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the leading edge radius of a 4-digit thickness
!   distribution.

  REAL,INTENT(IN):: toc       ! maximum value of t/c

  REAL:: rle
  REAL,PARAMETER:: A=1.1019  ! Eq. 6.3, p.114 in Abbott & von Doenhoff
!-------------------------------------------------------------------------------
  rle=A*toc**2
  RETURN
END Function LeadingEdgeRadius4    ! -------------------------------------------

!+
PURE FUNCTION LeadingEdgeRadius4M(toc,leIndex) RESULT(rle)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the leading edge radius of a 4-digit-modified thickness
!   distribution from the thickness and the le index.
! NOTE - leIndex is coded as REAL, although it is usually integer.

  REAL,INTENT(IN):: toc       ! maximum value of t/c
  REAL,INTENT(IN):: leIndex   ! leading edge index

  REAL:: rle
  REAL,PARAMETER:: A=1.1019/36.0  ! see p.117 in Abbott & von Doenhoff
!-------------------------------------------------------------------------------
  rle=A*(toc*leIndex)**2
  RETURN
END Function LeadingEdgeRadius4M    ! ------------------------------------------

!+
FUNCTION LeadingEdgeRadius6(family,toc) RESULT(rle)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the leading edge radius of a 6- or 6A-series thickness
!   distribution.
USE SplineProcedures,ONLY: PClookup,FMMspline

  INTEGER,INTENT(IN):: family
  REAL,INTENT(IN):: toc       ! maximum value of t/c
  REAL:: rle
  
  REAL,DIMENSION(201):: xt,yt
  REAL,DIMENSION(401):: sLocal,xLocal,yLocal,xpLocal,ypLocal
  REAL:: xpp, yp
!
!  REAL,PARAMETER,DIMENSION(7,8):: R = RESHAPE( (/                       &
!    0.0,0.297,0.631,1.087,1.594,2.120,2.650,                            &
!    0.0,0.256,0.579,1.040,1.590,2.208,2.884,                            &
!    0.0,0.240,0.552,1.000,1.505,1.960,2.500,                            &
!    0.0,0.223,0.530,0.952,1.435,1.955,2.550,                            &
!    0.0,0.223,0.530,0.952,1.435,1.955,2.550,                            &
!    0.0,0.265,0.608,1.071,1.630,2.120,2.650,                            &
!    0.0,0.246,0.579,0.994,1.561,2.208,2.884,                            &
!    0.0,0.229,0.524,0.922,1.446,1.955,2.550 /), (/7,8/) )
!  REAL,PARAMETER,DIMENSION(7):: TC = (/0.0,0.06,0.09,0.12,0.15,0.18,0.21 /)  
!-------------------------------------------------------------------------------
  CALL SetSixDigitPoints(family,toc, xt,yt)
  CALL ParametrizeAirfoil(xt,yt, xt,-yt, sLocal,xLocal,yLocal)
  CALL FMMSpline(sLocal,xLocal,xpLocal)
  CALL FMMSpline(sLocal,yLocal,ypLocal)
  CALL PClookup(slocal,xLocal,xpLocal, sLocal(201), FPP=xpp)
  CALL PClookup(slocal,yLocal,ypLocal, sLocal(201), FP=yp)
  rle=yp*yp/xpp
  
  RETURN
END Function LeadingEdgeRadius6    ! -------------------------------------------

!+
SUBROUTINE MeanLine2(cmax,xmaxc,x,ym,ymp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute displacement and slope of the NACA 2-digit mean line.
! NOTE - E.g., a NACA Mean Line 64 has cmax=0.06 and xmax=0.4
! REF - Eq 6.4, p.114 of Abbott and von Doenhoff

IMPLICIT NONE

  REAL,INTENT(IN):: cmax   ! max. camber as fraction of chord
  REAL,INTENT(IN):: xmaxc  ! fraction of chord where the camber is max.
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(OUT),DIMENSION(:):: ym
  REAL,INTENT(OUT),DIMENSION(:):: ymp

  INTEGER:: k,n
  REAL:: slope,slope1,slope2
  REAL:: theta
!-------------------------------------------------------------------------------
  slope1=2.0*cmax/xmaxc
  slope2=-2.0*cmax/(1.0-xmaxc)
  n=SIZE(x)
  DO k=1,n
    IF (x(k) < xmaxc) THEN
      theta=x(k)/xmaxc
      slope=slope1
    ELSE
      theta=(1.0-x(k))/(1.0-xmaxc)
      slope=slope2
    END IF
    ym(k)=theta*(2.0-theta)
    ymp(k)=slope*(1.0-theta)
  END DO

  ym(1:n)=cmax*ym(1:n)
  RETURN 
END Subroutine MeanLine2   ! ---------------------------------------------------

!+
SUBROUTINE MeanLine3(cl,xmaxc,x,ym,ymp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute displacement and slope of the NACA 3-digit mean line
! EXAMPLE - A 210 mean line has design CL=0.3 and max camber at 5 percent
!   (first digit is 2/3 CL in tenths; 
!    second digit is twice x of max camber in percent chord;
!    third digit, zero, indicates non-reflexed )
! REF - Eq. 6.6, p.115 of Abbott and von Doenhoff
! NOTE - The constant 1.8 appears in the final calculation. This comes
!   from the fact that the basic equation is given for cl=0.3 and that
!   there is a 6 in the denominator. So, one finally multiplies by
!   (cl/0.3)/6

  REAL,INTENT(IN):: cl ! design lift coefficient  
  REAL,INTENT(IN):: xmaxc  ! x-coor of maximum camber
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(OUT),DIMENSION(:):: ym
  REAL,INTENT(OUT),DIMENSION(:):: ymp

  INTEGER:: k,n
  REAL:: k1  ! related to y at x=r
  REAL:: r  ! fraction of chord where 2nd deriv goes to zero 
  REAL:: xx
!-------------------------------------------------------------------------------
  n=SIZE(x)
  CALL GetRk1(xmaxc,r,k1)   ! computes r and k1
!  WRITE(*,*) "3-digit, r,k1=", r,k1

  DO k=1,n
    xx=x(k)
    IF (xx < r) THEN
      ym(k)=xx*(xx*(xx-3.0*r)+r*r*(3.0-r))
      ymp(k)=3.0*xx*(xx-r-r)+r*r*(3.0-r)
    ELSE
      ym(k)=r*r*r*(1.0-xx)  
      ymp(k)=-r*r*r
    END IF
  END DO

  ym(1:n)=(k1*cl/1.8)*ym(1:n)
  ymp(1:n)=(k1*cl/1.8)*ymp(1:n)

  RETURN
END Subroutine MeanLine3   ! ---------------------------------------------------

!+
SUBROUTINE MeanLine3Reflex(cl,xmaxc, x,ym,ymp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute displacement and slope of the NACA 
!   3-digit reflex mean line
! NOTE - The constant 1.8 appears in the final calculation. This comes
!   from the fact that the basic equation is given for cl=0.3 and that
!   there is a 6 in the denominator. So, one finally multiplies by
!   (cl/0.3)/6
! REF - p.8 of NASA Technical Memorandum 4741 and NACA Report 537.

  REAL,INTENT(IN):: cl  ! design lift coefficient
  REAL,INTENT(IN):: xmaxc ! x-coor of maximum camber
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(OUT),DIMENSION(:):: ym
  REAL,INTENT(OUT),DIMENSION(:):: ymp

  INTEGER:: k,n
  REAL:: k1,k21
  REAL:: mr3   ! (1.0-r)**3
  REAL:: r
  REAL:: r3
  REAL:: xx
!---------------------------------------------------------------------------
  n=SIZE(x)
  CALL GetRk1k2(xmaxc,r,k1,k21)   ! computes r,k1,k2
!  WRITE(*,*) "3-digit-reflex, r,k1,k21=", r,k1,k21
  
  r3=r**3
  mr3=(1.0-r)**3
 
  DO k=1,n
    xx=x(k)
    IF (xx < r) THEN
      ym(k)=(xx-r)**3-k21*mr3*xx-xx*r3+r3
      ymp(k)=3.0*(xx-r)**2-k21*mr3-r3
    ELSE
      ym(k)=k21*(xx-r)**3-k21*mr3*xx-xx*r3+r3
      ymp(k)=3.0*k21*(xx-r)**2-k21*mr3-r3
    END IF
  END DO

  ym(1:n)=(k1*cl/1.8)*ym(1:n)
  ymp(1:n)=(k1*cl/1.8)*ymp(1:n)

  RETURN
END Subroutine MeanLine3Reflex   ! ---------------------------------------------

!+
SUBROUTINE MeanLine6(a,cl,x,ym,ymp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the displacement and slope of a 6-series mean line.
! REF - Eq. 4-26 and 4-27, p.74 of Abbott and von Doenhoff.

  REAL,INTENT(IN):: a  ! chordwise extent of uniform loading 0<=a<=1
  REAL,INTENT(IN):: cl ! design CL of the mean line
  REAL,INTENT(IN), DIMENSION(:):: x   ! fraction of chord
  REAL,INTENT(OUT),DIMENSION(:):: ym   ! y-coor of meanline
  REAL,INTENT(OUT),DIMENSION(:):: ymp   ! dy/dx of mean line

  REAL,PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, HALF=0.5, FOURTH=0.25
  REAL,PARAMETER:: PI=3.141592654, TWOPI=TWO*PI, EPS= 1E-7

  REAL:: amx
  REAL:: g,h
  INTEGER:: k,n
  REAL:: oma
  REAL:: omx
  REAL:: term1,term2, term1p,term2p
  REAL:: xx
!-------------------------------------------------------------------------------
  n=SIZE(x)
  ym(1:n) = ZERO
  ymp(1:n)= ZERO

  oma=ONE-a
  IF (ABS(oma) < EPS) THEN    
    DO k=1,n                    ! a==1
      xx=x(k)
      omx= ONE-xx
      IF (xx<EPS .OR. omx<EPS) THEN
        ym(k)=ZERO
        ymp(k)=ZERO
        CYCLE
      END IF
      ym(k)=omx*LOG(omx)+xx*LOG(xx)
      ymp(k)=LOG(omx)-LOG(xx)
    END DO
    ym(1:n)=-ym(1:n)*cl*(FOURTH/PI)
    ymp(1:n)=ymp(1:n)*cl*(FOURTH/PI)
    RETURN
  END IF

  DO k=1,n 
    xx=x(k)
    omx= ONE-xx
    IF (xx<EPS .OR. ABS(omx)<EPS) THEN
      ym(k)=ZERO
      ymp(k)=ZERO
      CYCLE
    END IF

    IF (ABS(a) < EPS) THEN                               ! a=0
      g=-FOURTH
      h=-HALF
    ELSE
      g=-(a*a*(HALF*LOG(a)-FOURTH)+FOURTH)/oma
      h=g+(HALF*oma*oma*LOG(oma)-FOURTH*oma*oma)/oma
    END IF

    amx= a-xx
    IF (ABS(amx) < EPS) THEN
      term1=ZERO
      term1p=ZERO
    ELSE
      term1=amx*amx*(TWO*LOG(ABS(amx))-ONE)
      term1p=-amx*LOG(ABS(amx))
    END IF

    term2=omx*omx*(ONE-TWO*LOG(omx))
    term2p=omx*LOG(omx)

    ym(k)=FOURTH*(term1+term2)/oma - xx*LOG(xx) + g - h*xx
    ymp(k)=(term1p+term2p)/oma - ONE - LOG(xx) - h
  END DO

  ym(1:n)=cl*ym(1:n)/(TWOPI*(a+ONE))
  ymp(1:n)=cl*ymp(1:n)/(TWOPI*(a+ONE))

  RETURN
END Subroutine MeanLine6  ! ----------------------------------------------------

!+
SUBROUTINE MeanLine6M(cl,x,ym,ymp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the displacement and slope of a modified six-series mean 
!   line as designed for use with the 6A-series profiles.
! REF - p.210 of NACA Report 903 (Larry Loftin).
! NOTE - a is not an input variable; it is always 0.8 for this mean line.

IMPLICIT NONE

  REAL,INTENT(IN):: cl   ! the design lift coefficient
 
  REAL,INTENT(IN),DIMENSION(:):: x   ! fraction of chord
  REAL,INTENT(OUT),DIMENSION(:):: ym   ! y-coor of meanline
  REAL,INTENT(OUT),DIMENSION(:):: ymp   ! dy/dx of mean line

  REAL,PARAMETER:: A = 0.8   ! chordwise extent of uniform loading 0<=a<=1
  REAL,PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, HALF=0.5, FOURTH=0.25
  REAL,PARAMETER:: PI=3.141592654, TWOPI=TWO*PI, EPS= 1E-7

  REAL:: amx
  REAL:: g,h
  INTEGER:: k,n
  REAL:: oma
  REAL:: omx
  REAL:: term1,term2, term1p,term2p
  REAL:: teSlope
  REAL:: xx
!-------------------------------------------------------------------------------
  n=SIZE(x)
  ym(1:n)=ZERO
  ymp(1:n)=ZERO

  oma=ONE-A
  teSlope=-0.24521*cl

  DO k=1,n
    xx=x(k)

    IF (xx<EPS) THEN
      ym(k)=ZERO
      ymp(k)=ZERO
      CYCLE
    END IF

    omx= ONE-xx
    IF (omx < EPS) THEN
      ym(k)=ZERO
      ymp(k)=teSlope
      CYCLE
    END IF

    g=-(A*A*(HALF*LOG(A)-FOURTH)+FOURTH)/oma
    h=g+(HALF*oma*oma*LOG(oma)-FOURTH*oma*oma)/oma

    amx= A-xx
    IF (ABS(amx) < EPS) THEN
      term1=ZERO
      term1p=ZERO
    ELSE
      term1=amx*amx*(TWO*LOG(ABS(amx))-ONE)
      term1p=-amx*LOG(ABS(amx))
    END IF

    term2=omx*omx*(ONE-TWO*LOG(omx))
    term2p=omx*LOG(omx)

    ym(k)=FOURTH*(term1+term2)/oma - xx*LOG(xx) + g - h*xx
    ymp(k)=(term1p+term2p)/oma - ONE - LOG(xx) - h
  END DO

  ym(1:n)=ym(1:n)*cl*0.97948/(TWOPI*(A+ONE))
  ymp(1:n)=ymp(1:n)*cl*0.97948/(TWOPI*(A+ONE))

  DO k=1,n   ! this is the special treatment...
    IF (x(k) > 0.86) THEN      
      ym(k)=teSlope*(x(k)-1.0)
      ymp(k)=teSlope
    END IF
  END DO

  RETURN
END Subroutine MeanLine6M  ! ---------------------------------------------------

!+
SUBROUTINE ParametrizeAirfoil(xupper,yupper,xlower,ylower, s,x,y)
! ------------------------------------------------------------------------------
! PURPOSE - Define the shape of an airfoil parametrically using s, the
!   inscribed arc length, s. s=0 at the trailing edge, then increases by
!   moving forward along the upper surface, then around the nose and rearward
!   along the lower surface. The final value of s is somewhat greater than 2.
  REAL,INTENT(IN),DIMENSION(:):: xupper,yupper,xlower,ylower
  REAL,INTENT(OUT),DIMENSION(:):: s,x,y

  INTEGER:: k,nn
  INTEGER:: nupper,nlower   ! don't have to be the same
!-------------------------------------------------------------------------------
  nupper=SIZE(xupper)
  nlower=SIZE(xlower)

  IF (SIZE(yupper)<nupper .OR. SIZE(ylower)<nlower) THEN
    WRITE(*,*) "DISASTER #1"
    STOP
  END IF

  nn=nupper+nlower-1
  IF (SIZE(s)<nn .OR. SIZE(x)<nn .OR. SIZE(y)<nn) THEN
    WRITE(*,*) "DISASTER #2"
    STOP
  END IF

  x(1:nupper)=xupper(nupper:1:-1)
  y(1:nupper)=yupper(nupper:1:-1)
  x(nupper+1:nn)=xlower(2:nlower)
  y(nupper+1:nn)=ylower(2:nlower)

  s(1)=0.0
  DO k=2,nn
    s(k)=s(k-1) + SQRT((x(k)-x(k-1))**2 + (y(k)-y(k-1))**2)
  END DO

  RETURN
END Subroutine ParametrizeAirfoil   ! ------------------------------------------

!+
PURE FUNCTION Polynomial(c,x) RESULT(f)
! ------------------------------------------------------------------------------
! PURPOSE - Evaluate a polynomial with real coefficients at x
!   f= c(1) + c(2)*x + c(3)*x**2 + c(4)*x**3 + c(5)*x**4 + ...
!   Straightforward application of Horner's rule.

  REAL,INTENT(IN),DIMENSION(:):: c
  REAL,INTENT(IN):: x
  REAL:: f

  REAL:: ff

  INTEGER:: j,n
!-------------------------------------------------------------------------------
  n=SIZE(c)
  ff=c(n)
  DO j=n-1,1,-1
    ff = ff*x + c(j)
  END DO

  f=ff
  RETURN
END Function Polynomial   ! ----------------------------------------------------

!+
PURE FUNCTION ScaleFactor(family,tc) RESULT(s)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the scale factor used to multiply the epsilon and psi
!   functions to get the proper thickness/chord for the mapping of a 
!   six-series airfoil (or 6A).
! REF - Discussed in AIAA-2001-5235.

  INTEGER,INTENT(IN):: family  ! =1 63; =2 64; =3 65; =4 66; =5 67;
                               ! =6 63A; =7 64A; =8 65A 
  REAL,INTENT(IN):: tc   ! desired t/c (fraction, not percent)
  REAL:: s   ! the scale factor to use

!... a column for each of the 6- and 6a-series families
  REAL,PARAMETER,DIMENSION(5,8):: COEFF = RESHAPE( (/      &
    0.0, 8.1827699, 1.3776209,  -0.092851684, 7.5942563,   &
    0.0, 4.6535511, 1.038063,   -1.5041794,   4.7882784,   &
    0.0, 6.5718716, 0.49376292,  0.7319794,   1.9491474,   &
    0.0, 6.7581414, 0.19253769,  0.81282621,  0.85202897,  &
    0.0, 6.627289,  0.098965859, 0.96759774,  0.90537584,  &
    0.0, 8.1845925, 1.0492569,   1.31150930,  4.4515579,   &
    0.0, 8.2125018, 0.76855961,  1.4922345,   3.6130133,   &
    0.0, 8.2514822, 0.46569361,  1.50113018,  2.0908904 /),  (/5,8/) )
!-------------------------------------------------------------------------------
  IF (family <1 .OR. family>8 .OR. tc <=0.0) THEN
    s=0.0
    RETURN
  END IF
  s=Polynomial(COEFF(:,family), tc)
  RETURN
END Function ScaleFactor   ! ---------------------------------------------------

!+
PURE SUBROUTINE SetSixDigitPoints(family,tc,xt,yt)
! ------------------------------------------------------------------------------
! PURPOSE - Set the data points that define a 6-series or 6A-series
!   thickness distribution. User has no control of the spacing of xt - it
!   depends on phi,eps,psi

USE EpsilonPsi
  INTEGER,INTENT(IN):: family   ! =1 63; =2 64; ... ; =6 63A; =7 64A; =8 65A
  REAL,INTENT(IN):: tc    ! max value of t/c
  REAL,INTENT(OUT),DIMENSION(:):: xt
  REAL,INTENT(OUT),DIMENSION(:):: yt

  REAL,PARAMETER:: a=1.0
  REAL,DIMENSION(201):: phi,eps,psi
  INTEGER:: i
  INTEGER,PARAMETER:: NP=201
  REAL:: sf
  COMPLEX,DIMENSION(201):: z,zprime,zeta,zfinal
!-------------------------------------------------------------------------------
  phi= (/ (REAL(i),i=0,NP-1)   /)
  phi=(3.14159265/REAL(NP-1))*phi
  SELECT CASE(family)   ! copy the appropriate epsilon and psi functions from
                        ! the data sets in module EpsilonPsi
    CASE(1)
      eps=EPS1
      psi=PSI1
    CASE(2)
      eps=EPS2
      psi=PSI2
    CASE(3)
      eps=EPS3
      psi=PSI3
    CASE(4)
      eps=EPS4
      psi=PSI4
    CASE(5)
      eps=EPS5
      psi=PSI5
    CASE(6)
      eps=EPS6
      psi=PSI6
    CASE(7)
      eps=EPS7
      psi=PSI7
    CASE(8)
      eps=EPS8
      psi=PSI8
!    CASE DEFAULT
      
  END SELECT

  sf=ScaleFactor(family,tc)
  eps=sf*eps
  psi=sf*psi

!... Now use this scaled set of epsilon and psi functions to perform the
!.   conformal mapping of the circle z into the scaled airfoil zfinal.
!.   Return the real and imaginary parts as xt and yt.
  z=a*EXP(CMPLX(psi(1),phi))
  zprime=z*EXP(CMPLX(psi-psi(1),-eps))
  zeta=zprime+a*a/zprime
  zfinal=(zeta(1)-zeta)/ABS(zeta(SIZE(zeta))-zeta(1))
  xt=REAL(zfinal)
  yt=-AIMAG(zfinal)

  RETURN
END Subroutine SetSixDigitPoints   ! -------------------------------------------

!+
SUBROUTINE Thickness4(toc,x,y,yp,ypp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the values of y, dy/dx, d2y/dx2 at each point of x

IMPLICIT NONE
  REAL,INTENT(IN):: toc
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(OUT),DIMENSION(:):: y
  REAL,INTENT(OUT),DIMENSION(:),OPTIONAL:: yp
  REAL,INTENT(OUT),DIMENSION(:),OPTIONAL:: ypp

  REAL,PARAMETER:: A0 =  0.2969
  REAL,PARAMETER:: A1 = -0.1260
  REAL,PARAMETER:: A2 = -0.3516, A22=A2+A2
  REAL,PARAMETER:: A3 =  0.2843, A33=A3+A3+A3
  REAL,PARAMETER:: A4 = -0.1015, A42=A4+A4, A44=A42+A42
  INTEGER:: k,n
  REAL:: srx,xx
!-------------------------------------------------------------------------------
  n=SIZE(x)
  IF (Present(yp))  yp(1)=1E22
  IF (Present(ypp)) ypp(1)=-1E22
  DO k=1,n    ! based on t/c=0.2
    xx=x(k)
    IF (xx==0.0) THEN
      y(k)=0.0
      IF (Present(yp))   yp(k)=1E22
      IF (Present(ypp)) ypp(k)=1E22
    ELSE
      srx=SQRT(xx)
      y(k)=A0*srx+xx*(A1+xx*(A2+xx*(A3+xx*A4)))
      IF (Present(yp)) &
        yp(k)=0.5*A0/srx + A1 + xx*(A22+xx*(A33+xx*A44))
      IF (Present(ypp)) &
        ypp(k)=-0.25*A0/(srx*srx*srx) + A22 + 6.0*xx*(A3+xx*A42)
    END IF
  END DO

  y(1:n)=(5.0*toc)*y(1:n)       ! convert to correct t/c
  IF (Present(yp)) yp(1:n)=(5.0*toc)*yp(1:n)
  IF (Present(ypp)) ypp(1:n)=(5.0*toc)*ypp(1:n)

  RETURN
END Subroutine Thickness4   ! --------------------------------------------------

!+
SUBROUTINE Thickness4M(toc,leIndex,xmaxt,x, y,yp,ypp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the thickness distribution (and derivatives) of a
!    NACA 4-digit-modified section.
! NOTE - first digit after dash is l.e. index; 2nd is loc of max thickness

IMPLICIT NONE
  REAL,INTENT(IN):: toc
  REAL,INTENT(IN):: leIndex  ! leading edge index
  REAL,INTENT(IN):: xmaxt  ! x-coor of maximum thickness (fraction of chord)
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(OUT),DIMENSION(:):: y
  REAL,INTENT(OUT),DIMENSION(:),OPTIONAL:: yp,ypp

  REAL:: a0,a1,a2,a3
  REAL:: a22,a33
  REAL:: d0,d1,d2,d3
  REAL:: d22,d33
  INTEGER:: k,n
  REAL:: omxmt,omxmsq
  REAL:: rle
  REAL:: srx
  REAL:: xx
!-------------------------------------------------------------------------------
  n=SIZE(x)

  d1=CalculateD1(xmaxt)
  rle=LeadingEdgeRadius4M(toc,leIndex)

  a0=(0.2/toc)*SQRT(rle+rle)
  a3= Cala3(a0, d1, xmaxt)
  a2= Cala2(a0, a3, xmaxt)
  a1= Cala1(a0, a2, a3, xmaxt)
  a22=a2+a2
  a33=a3+a3+a3

  omxmt= 1.0 - xmaxt !  calculate the "d" constants
  omxmsq= omxmt * omxmt

  d3= ((3. * d1) - (.588 / omxmt)) / (3. * omxmsq)
  d2= (-1.5 * omxmt * d3) - ((.5 * d1) / omxmt)
  d0= .002
  d22=d2+d2
  d33=d3+d3+d3

!  WRITE(*,*) "l.e. radius=", rle       ! remove after checkout
!  WRITE(*,*) "A:", a0,a1,a2,a3         ! remove after checkout
!  WRITE(*,*) "D:", d0,d1,d2,d3         ! remove after checkout

  DO k=1,n   ! based on t/c=0.2
    xx=x(k)
    IF (xx==0.0) THEN
      y(k)=0.0
      IF (Present(yp))   yp(k)=1E20
      IF (Present(ypp)) ypp(k)=-1E20
      CYCLE
    END IF
    IF (x(k) < xmaxt) THEN
      srx= SQRT(xx)
      y(k)=A0*srx + xx*(A1 + xx*(A2 + xx*A3))
      IF (Present(yp))   yp(k)=0.5*A0/srx + A1 + xx*(A22+xx*A33)
      IF (Present(ypp)) ypp(k)=-0.25*A0/(srx*srx*srx) + A22 + A33*xx
    ELSE
      xx=1.0-x(k)
      y(k)=D0 + xx*(D1 + xx*(D2 + xx*D3))
      IF (Present(yp))   yp(k)=D1+xx*(D22+xx*D33)
      IF (Present(ypp)) ypp(k)=D22+2.0*xx*D33
    END IF
  END DO

  y(1:n)=5.0*toc*y(1:n)              ! convert to correct t/c
  IF (Present(yp))  yp(1:n)=5.0*toc*yp(1:n)
  IF (Present(ypp)) ypp(1:n)=5.0*toc*ypp(1:n)

  RETURN
END Subroutine Thickness4M   ! -------------------------------------------------

!+
SUBROUTINE Thickness4sharpTE(toc,x,y,yp,ypp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the values of y, dy/dx, d2y/dx2 at each point of x
!   The equation used gives an airfoil that is very similar to the four digit
!   airfoil profile, but has a sharp trailing edge.
!  This airfoil does not have any NACA blessing, so it is strictly unofficial.
!  Lots of folks want to compute a four-digit section, but are using a theory
!  that assumes a zero trailing edge.

IMPLICIT NONE
  REAL,INTENT(IN):: toc
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(OUT),DIMENSION(:):: y
  REAL,INTENT(OUT),DIMENSION(:),OPTIONAL:: yp
  REAL,INTENT(OUT),DIMENSION(:),OPTIONAL:: ypp

  REAL,PARAMETER:: A0 =  0.2969
  REAL,PARAMETER:: A1 = -0.1260
  REAL,PARAMETER:: A2 = -0.3516, A22=A2+A2
  REAL,PARAMETER:: A3 =  0.2843, A33=A3+A3+A3
  REAL,PARAMETER:: A4 = -0.1036, A42=A4+A4, A44=A42+A42   ! slightly changed
  INTEGER:: k,n
  REAL:: srx,xx
!-------------------------------------------------------------------------------
  n=SIZE(x)
  IF (Present(yp))  yp(1)=1E22
  IF (Present(ypp)) ypp(1)=-1E22
  DO k=1,n    ! based on t/c=0.2
    xx=x(k)
    IF (xx==0.0) THEN
      y(k)=0.0
      IF (Present(yp))   yp(k)=1E22
      IF (Present(ypp)) ypp(k)=1E22
    ELSE
      srx=SQRT(xx)
      y(k)=A0*srx+xx*(A1+xx*(A2+xx*(A3+xx*A4)))
      IF (Present(yp)) &
        yp(k)=0.5*A0/srx + A1 + xx*(A22+xx*(A33+xx*A44))
      IF (Present(ypp)) &
        ypp(k)=-0.25*A0/(srx*srx*srx) + A22 + 6.0*xx*(A3+xx*A42)
    END IF
  END DO

  y(1:n)=(5.0*toc)*y(1:n)       ! convert to correct t/c
  IF (Present(yp)) yp(1:n)=(5.0*toc)*yp(1:n)
  IF (Present(ypp)) ypp(1:n)=(5.0*toc)*ypp(1:n)

  RETURN
END Subroutine Thickness4sharpTE   ! -------------------------------------------

!+
SUBROUTINE Thickness6(family,toc,x, y,yp,ypp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the values of y, dy/dx, d2y/dx2 at each point of x
!   by interpolation in the tables computed by xxxx.

USE SplineProcedures,ONLY: FMMspline,PClookup,SplineZero
  INTEGER,INTENT(IN):: family
  REAL,INTENT(IN):: toc
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(OUT),DIMENSION(:):: y
  REAL,INTENT(OUT),DIMENSION(:),OPTIONAL:: yp
  REAL,INTENT(OUT),DIMENSION(:),OPTIONAL:: ypp

  INTEGER:: errCode
!  REAL:: errMax
  INTEGER:: k
  REAL:: sx
  REAL,PARAMETER:: TOL=1E-6
  REAL,DIMENSION(201):: xt,yt
  REAL,DIMENSION(401):: xLocal,yLocal,sLocal,  xpLocal,ypLocal
  REAL,DIMENSION(201):: sLower,xLower,yLower, xpLower,ypLower  !!,spLower
  REAL,DIMENSION(SIZE(x)):: dxds,d2xds2, dyds,d2yds2
!-------------------------------------------------------------------------------
  IF (Present(yp))   yp=0.0
  IF (Present(ypp)) ypp=0.0
  CALL SetSixDigitPoints(family,toc, xt,yt)

  CALL ParametrizeAirfoil(xt,yt, xt,-yt, sLocal,xLocal,yLocal)

! now fit splines to xLocal vs. sLocal  and  yLocal vs. slocal
  CALL FMMSpline(sLocal,xLocal,xpLocal)
  CALL FMMSpline(sLocal,yLocal,ypLocal)

  sLower=sLocal(201:401)
  xLower=xLocal(201:401)
  yLower=-yLocal(201:401)
  xpLower=xpLocal(201:401)
  ypLower=-ypLocal(201:401)
!  spLower(1)=0.0                      !  spLower is ds/dx
!  spLower(2:201)=1.0/xpLower(2:201)

  DO k=1,SIZE(x)
    CALL SplineZero(sLower,xLower,xpLower, x(k),TOL,sx,errCode)
    IF (errCode /= 0) WRITE(*,*) "Non-zero errCode in Thickness6"
!!!  CALL PClookup(xLower,sLower,spLower, x,sx)
    CALL PClookup(sLower,xLower,xpLower, sx, FP=dxds(k), FPP=d2xds2(k))
    CALL PClookup(slower,yLower,ypLower, sx, y(k), dyds(k), d2yds2(k))
  END DO
  
!  errMax=MAXVAL(ABS(xdum-x))
!  WRITE(*,*) "max error in x=", errMax
!  WRITE(IDBG,*) "max error in x=", errMax

  IF (Present(yp)) THEN    
    WHERE (dxds == 0.0)
      yp(1:SIZE(x))=0.0    
    ELSEWHERE
      yp(1:SIZE(x))=dyds(:)/dxds(:)
    END WHERE
  END IF

  IF (Present(ypp)) THEN
    WHERE (d2xds2 == 0.0)
      ypp(1:SIZE(x))=0.0
    ELSEWHERE
      ypp(1:SIZE(x))=d2yds2(:)/d2xds2(:)
    END WHERE
  END IF

  RETURN
END Subroutine Thickness6   ! --------------------------------------------------

!+

END Module NacaAuxiliary   ! ===================================================
