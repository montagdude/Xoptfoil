!+
MODULE SplineProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Collect several procedures associated with spline functions that
!  are used in the NACA airfoil calculations.
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
! DATE   PERSON      STATEMENT OF CHANGES
! 09Oct01  RLC   0.5   Original coding
! 24Nov01  RLC   0.6   Made a,fa... PRIVATE
! 05Dec01  RLC   0.7   Added material from Cubics and IntUtil
! 26Dec01  RLC   0.8   Made function values optional in PClookup and elsewhere
!  4Jan02  RLC   1.0   Final cleanup for release of PDAS 7
! 16Jan09  RLC   1.1   Additional cosmetic improvements


  IMPLICIT NONE
!-------------------------------------------------------------------------------
  CHARACTER(LEN=*),PARAMETER,PUBLIC:: &
    SPLINE_PROCEDURES_VERSION="1.1 (16 January 2009)"

  PRIVATE:: EvaluateCubic
  PUBLIC:: EvaluateCubicAndDerivs
  PUBLIC:: FMMspline
  PRIVATE:: InterpolatePolynomial
  PRIVATE:: LookUp
  PUBLIC:: PClookup
  PUBLIC:: SplineZero
  PUBLIC:: TableLookup
  PRIVATE:: Zeroin

!... Module-wide variables set in SplineZero and used by EvaluateCubic
  REAL,PRIVATE:: a,fa,fpa
  REAL,PRIVATE:: b,fb,fpb


CONTAINS

!+
PURE FUNCTION EvaluateCubic(u) RESULT(fu)
! ------------------------------------------------------------------------------
! PURPOSE - Evaluate a cubic polynomial defined by the function and the
!   1st derivative at two points
  REAL,INTENT(IN):: u   ! point where function is to be evaluated
  REAL:: fu                    ! computed value of f(u)

!  REAL:: a,fa,fpa   ! a, f(a), f'(a)  at first point
!  REAL:: b,fb,fpb   ! b, f(b), f'(b)  at second point
  

  REAL:: d,t,p
!-------------------------------------------------------------------------------
  d=(fb-fa)/(b-a)
  t=(u-a)/(b-a)
  p=1.0-t

  fu = p*fa + t*fb - p*t*(b-a)*(p*(d-fpa)-t*(d-fpb))
  RETURN
END Function EvaluateCubic   ! -------------------------------------------------

!+
PURE SUBROUTINE EvaluateCubicAndDerivs(a,fa,fpa, b,fb,fpb, u, f,fp,fpp,fppp)
! ------------------------------------------------------------------------------
! PURPOSE - Evaluate a cubic polynomial and its 1st, 2nd, and 3rd
!   derivatives at a specified point. The cubic is defined by the function
!   and the 1st derivative at two points. The cubic is mapped onto [0,1]
!   for some savings in computation. If no derivatives are wanted, it is
!   probably preferable to use functions EvaluateCubic above.

  REAL,INTENT(IN):: a,fa,fpa   ! a, f(a), f'(a)  at first point
  REAL,INTENT(IN):: b,fb,fpb   ! b, f(b), f'(b)  at second point
  REAL,INTENT(IN):: u   ! point where function is to be evaluated

  REAL,INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp   ! f(u),f'(u),f''(u),f'''(u)
  
  REAL,PARAMETER:: ONE=1.0, TWO=2.0, THREE=3.0, SIX=6.0
! the "magic" matrix. This is how you declare a constant 4x4 matrix
  REAL,PARAMETER,DIMENSION(4,4):: MAGIC = RESHAPE(                      &
     & (/2.0, -3.0,  0.0,  1.0, -2.0,  3.0,  0.0,  0.0,                 &
     &   1.0, -2.0,  1.0,  0.0,  1.0, -1.0,  0.0,  0.0/), (/4,4/) )
  REAL,DIMENSION(4)::coef,rhs
  REAL:: h,t
! ------------------------------------------------------------------------------
  rhs(1)=fa
  rhs(2)=fb
  rhs(3)=fpa*(b-a)
  rhs(4)=fpb*(b-a)
  coef=MATMUL(MAGIC,rhs)

! CAUTION - these are not the coefficients of the cubic in the original 
! coordinates. This is the cubic on [0,1] from the mapping t=(x-a)/(b-a). 
! That is why the h terms appear in the derivatives.       

  h=ONE/(b-a)
  t=(u-a)*h
  IF (Present(f)) THEN
    f=             coef(4) +     t*(coef(3) + t*(coef(2)   + t*coef(1)))
  END IF
  IF (Present(fp)) THEN
    fp=         h*(coef(3) + t*(TWO*coef(2) + t*THREE*coef(1)))
  END IF
  IF (Present(fpp)) THEN
    fpp=  h*h*(TWO*coef(2) +  t*SIX*coef(1))
  END IF
  IF (Present(fppp)) THEN
    fppp=h*h*h*SIX*coef(1)
  END IF
  RETURN
END Subroutine EvaluateCubicAndDerivs   ! -------------------------------------

!+
PURE SUBROUTINE FMMspline(x,y, yp)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the cubic spline with endpoint conditions chosen
!    by FMM  (from the book by Forsythe,Malcolm & Moler)
  REAL,INTENT(IN),DIMENSION(:):: x,y
  REAL,INTENT(OUT),DIMENSION(:):: yp

  INTEGER:: i,n
  REAL,ALLOCATABLE,DIMENSION(:):: dx,dy,delta,dd,alpha,beta,sigma
  REAL:: deriv1,deriv2
!-------------------------------------------------------------------------------
  n=SIZE(x)   ! y and yp must be at least this size
  IF (n<2) THEN
    yp(1)=0.0
    RETURN
  END IF
  ALLOCATE(dx(n-1), dy(n-1), delta(n-1))
  dx(:)=x(2:n)-x(1:n-1)
  dy(:)=y(2:n)-y(1:n-1)
  delta(:)=dy(:)/dx(:)
  IF (n==2) THEN
    yp(1)=delta(1)
    yp(2)=delta(1)
    DEALLOCATE(dx,dy,delta)
    RETURN
  END IF
  ALLOCATE(dd(n-2))
  dd(:)=delta(2:n-1)-delta(1:n-2)
  IF (n==3) THEN
    deriv2=dd(1)/(x(3)-x(1))
    deriv1=delta(1)-deriv2*dx(1)
    yp(1)=deriv1
    yp(2)=deriv1+deriv2*dx(1)
    yp(3)=deriv1+deriv2*(x(3)-x(1))
    DEALLOCATE(dx,dy,delta,dd)
    RETURN
  END IF
! This gets rid of the trivial cases n=1,2,3. Assume from here on n>3

  ALLOCATE(alpha(n),beta(n),sigma(n))
  alpha(1)=-dx(1)
  alpha(2:n-1)=2.0*(dx(1:n-2)+dx(2:n-1))
  DO i=2,n-1                                                    ! serial loop
    alpha(i)=alpha(i)-dx(i-1)*dx(i-1)/alpha(i-1)            ! fwd elimination
  END DO
  alpha(n)=-dx(n-1)-dx(n-1)*dx(n-1)/alpha(n-1)

  beta(1)=dd(2)/(x(4)-x(2)) - dd(1)/(x(3)-x(1))
  beta(1)=beta(1)*dx(1)*dx(1)/(x(4)-x(1))

  beta(2:n-1)=dd(1:n-2)

  beta(n)=dd(n-2)/(x(n)-x(n-2)) - dd(n-3)/(x(n-1)-x(n-3))
  beta(n)=-beta(n)*dx(n-1)*dx(n-1)/(x(n)-x(n-3)) 

  DO i=2,n                                                    ! serial loop
    beta(i)=beta(i)-dx(i-1)*beta(i-1)/alpha(i-1)          ! fwd elimination
  END DO

  sigma(n)=beta(n)/alpha(n)
  DO i=n-1,1,-1                                 ! reverse order serial loop
    sigma(i)=(beta(i)-dx(i)*sigma(i+1))/alpha(i)        ! back substitution
  END DO

  yp(1:n-1)=delta-dx*(sigma(1:n-1)+sigma(1:n-1)+sigma(2:n))
  yp(n)=yp(n-1)+dx(n-1)*3.0*(sigma(n)+sigma(n-1))

  DEALLOCATE(dx,dy,delta, alpha,beta,sigma, dd)

  RETURN
END SUBROUTINE FMMspline   ! ---------------------------------------------------

!+
PURE FUNCTION InterpolatePolynomial(x,y,u) RESULT(sum)
! ------------------------------------------------------------------------------
! PURPOSE -Compute the value of the interpolating polynomial thru
!   x- and y-arrays at the x-value of u, using Lagrange's equation.

  REAL,INTENT(IN),DIMENSION(:):: x,y   ! tables of coordinates
  REAL,INTENT(IN):: u   ! value of x-coordinate for interpolation
  REAL:: sum
  INTEGER:: i,j
  REAL:: fact
  REAL,DIMENSION(SIZE(x)):: du
!-------------------------------------------------------------------------------
  du(:)=u-x(:)
  sum=0.0
  DO j=1,SIZE(x)
    fact=1.0
    DO i=1,SIZE(x)
      IF (i /= j) fact=fact*du(i)/(x(j)-x(i))
    END DO
    sum=sum+y(j)*fact
  END DO
  RETURN
END Function InterpolatePolynomial   ! -----------------------------------------

!+
PURE FUNCTION Lookup(xtab,x) RESULT (i)
! ------------------------------------------------------------------------------
! PURPOSE - Search a sorted (increasing) array to find the interval
!  bounding a given number. If n is the size of the array a,
!  return 0 if number x is < a(1)
!  return n if x > a(n).
!  return i if  a(i) <= x < a(i+1).
!  If x is exactly equal to a(n), return n-1.

  REAL,INTENT(IN),DIMENSION(:)::  xtab    ! input array                        
  REAL,INTENT(IN):: x  ! input number              

  INTEGER:: i  ! index of interval such that xtab(i) <= x < xtab(i+1)
               ! =0 if x < xtab(1) and =n if x > xtab(i)
  INTEGER:: j,k,n
!-------------------------------------------------------------------------------
  n=SIZE(xtab)
  IF (n <= 0) THEN
    i=-1
    RETURN
  END IF

  IF (x < xtab(1)) THEN
    i=0
    RETURN
  END IF

  IF (x > xtab(n)) THEN
    i=n
    RETURN
  END IF

  i=1 
  j=SIZE(xtab)
  DO
    IF (j <= i+1) EXIT
    k=(i+j)/2                                     ! integer division
    IF (x < xtab(k)) THEN
      j=k
    ELSE
      i=k
    END IF      
  END DO

  RETURN
END Function Lookup   ! --------------------------------------------------------

!+
SUBROUTINE PClookup(x,y,yp, u, f,fp,fpp,fppp)
! ------------------------------------------------------------------------------
! PURPOSE - Interpolate in a cubic spline at one point 
  REAL,INTENT(IN),DIMENSION(:):: x,y,yp   ! defines the cubic spline
  REAL,INTENT(IN):: u   ! point where spline is to be evaluated
  REAL,INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp   ! f(u),f'(u), f''(u), f'''(u)

  REAL:: ud, a,fa,fpa, b,fb,fpb
  REAL:: z,zp,zpp,zppp
  INTEGER:: k
!-------------------------------------------------------------------------------
  k=LookUp(x,u)
  k=MAX(1, MIN(SIZE(x)-1,k) )
  a=x(k)
  fa=y(k)
  fpa=yp(k)
  b=x(k+1)
  fb=y(k+1)
  fpb=yp(k+1)
  ud=u
  CALL EvaluateCubicAndDerivs(a,fa,fpa, b,fb,fpb, ud, z,zp,zpp,zppp)
  IF(Present(f)) f=z
  IF(Present(fp)) fp=zp
  IF(Present(fpp)) fpp=zpp
  IF(Present(fppp)) fppp=zppp

  RETURN
END Subroutine PClookup   ! ----------------------------------------------------

!+
SUBROUTINE SplineZero(x,f,fp, fbar,tol, xbar,errCode) 
! ------------------------------------------------------------------------------
! PURPOSE - Find a value of x corresponding to a value of fbar of the cubic
!  spline defined by arrays x,f,fp. f is the value of the spline at x and fp
!  is the first derivative.
! NOTES - The spline is searched for an interval that crosses the specified
!   value. Then Brent's method is used for finding the zero.

  IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: x,f,fp
  REAL,INTENT(IN):: fbar
  REAL,INTENT(IN):: tol
  REAL,INTENT(OUT):: xbar
  INTEGER,INTENT(OUT):: errCode

  INTEGER:: k,n
  REAL,ALLOCATABLE,DIMENSION(:):: fLocal
  REAL,PARAMETER:: ZERO=0.0
!-------------------------------------------------------------------------------
  n=SIZE(x)

  DO k=1,n   ! look for an exact match. Could happen...
    IF (ABS(f(k)-fbar) < tol) THEN
      xbar=x(k)
      errCode=0
      RETURN
    END IF
  END DO

  ALLOCATE(fLocal(n))
  fLocal(1:n)=f(1:n)-fbar   ! look for a zero of fLocal

  DO k=2,n
   IF ( fLocal(k-1)*fLocal(k) < ZERO) EXIT
  END DO
  IF (k==n+1) THEN   ! no crossing could be found
    errCode=1
    DEALLOCATE(fLocal)
    RETURN
  END IF

  errCode=0
  a=x(k-1)    ! set the global variables for EvaluateCubic
  fa=fLocal(k-1)
  fpa=fp(k-1)
  b=x(k)
  fb=fLocal(k)
  fpb=fp(k)
  DEALLOCATE(fLocal)

  xbar=Zeroin(a,b,EvaluateCubic,tol)

  RETURN
END Subroutine SplineZero   ! --------------------------------------------------

!+
PURE FUNCTION TableLookup(x,y,order,u) RESULT(fu)
! ------------------------------------------------------------------------------
! PURPOSE - Use polynomial evaluation for table lookup. Find points ahead
!   and behind evaluation point to match order of interpolation desired.

  REAL,INTENT(IN),DIMENSION(:):: x,y   ! data tables
  INTEGER,INTENT(IN):: order   ! order of interpolation
                               ! (1=linear, 2=quadratic...)
  REAL,INTENT(IN):: u          ! x-coor where function is to be evaluated
  REAL:: fu                    ! interpolated function value
  INTEGER:: j,m
!-------------------------------------------------------------------------------
  m=MIN(order+1, SIZE(x))  ! number of points used for interpolating poly
  j=Lookup(x,u)
  j=j-(m/2-1)
  j=MIN(1+SIZE(x)-m, j)                       ! j+m-1 must not exceed SIZE(x)
  j=MAX(1,j)                                  ! j must be positive
                       ! use points j thru j+m-1 for interpolation (m points)
  fu=InterpolatePolynomial(x(j:j+m-1),y(j:j+m-1),u)
  RETURN
END Function TableLookup   ! ---------------------------------------------------

!+
      FUNCTION Zeroin (ax,bx,f,tol) RESULT(z)
! ------------------------------------------------------------------------------
! PURPOSE - Compute a zero of f in the interval (ax,bx)
! AUTHORS - Van Winjngaarden,Decker,Brendt, et al
!           Forsythe, Malcolm, & Moler, 
!           Press, et.al, Numerical Recipes
!           Ralph L. Carmichael, Public Domain Aeronautical Software

!     REVISIONS - DATE   PERSON      STATEMENT OF CHANGES


      IMPLICIT NONE
      REAL,INTENT(IN):: ax,bx   ! left and right endpoints of interval
      REAL,INTENT(IN):: tol     ! desired interval of uncertainity 
      
      REAL:: z

      INTERFACE
        FUNCTION F(x) RESULT(g)
          IMPLICIT NONE
          REAL,INTENT(IN):: x
          REAL:: g
        END Function F
      END INTERFACE    

      INTEGER,PARAMETER:: MAX_ITER=500
      INTEGER:: k
      REAL,PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0, HALF=0.5
      REAL:: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
!-------------------------------------------------------------------------------
      eps=EPSILON(ax)
!      tol1=ONE+eps

      a=ax   ! initialization
      b=bx
      fa=f(a)
      fb=f(b)   ! should test that fa and fb have opposite signs
      c=b
      fc=fb

      DO k=1,MAX_ITER   ! begin iteration

        IF ((fb > ZERO .AND. fc > ZERO) .OR. (fb < ZERO .AND. fc < ZERO)) THEN
          c=a
          fc=fa
          d=b-a
          e=d
        END IF

        IF (ABS(fc) < ABS(fb)) THEN
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        END IF
      
        tol1=TWO*eps*ABS(b) + HALF*tol   ! convergence test
        xm=0.5*(c-b)
        IF (ABS(xm)<=tol1 .OR. fb==0.0) THEN
          z=b
          RETURN   ! the only way out!
        END IF   

!.....is bisection necessary?...........................................
        IF (ABS(e)<=tol1 .AND. ABS(fa)>ABS(fb)) THEN
          s=fb/fa   ! is quadratic interpolation possible ?
          IF (a==c) THEN
            s=fb/fa   ! use linear interpolation
            p=TWO*xm*s
            q=ONE-s
          ELSE
            q=fa/fc   ! use inverse quadratic interpolation
            r=fb/fc
            s=fb/fa
            p=s*(TWO*xm*q*(q-r)-(b-a)*(r-ONE))
            q=(q-ONE)*(r-ONE)*(s-ONE)
          END IF

          IF (p > ZERO) q=-q   ! adjust signs
          p=ABS(p)

          IF (p+p < MIN(THREE*xm*q-ABS(tol1*q),ABS(e*q))) THEN
            e=d    ! use interpolation
            d=p/q
          ELSE
            d=xm   ! use bisection
            e=d
          END IF
        ELSE
          d=xm   ! use bisection
          e=d
        END IF


        a=b
        fa=fb
        IF (ABS(d)> tol1) THEN
          b=b+d
        ELSE
          b=b+SIGN(TOL1,XM)
        END IF
        fb=f(b)

      END DO

      z=b   ! but this is a bad return. Max iterations exceeded.
      RETURN
      END Function Zeroin   ! --------------------------------------------------

END Module SplineProcedures   ! ================================================
