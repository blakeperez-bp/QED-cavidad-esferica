SUBROUTINE gauleg(x1,x2,x,w,n)
INTEGER n
REAL*8 x1,x2,x(n),w(n)
DOUBLE PRECISION EPS        !  EPS is the relative precision.
PARAMETER (EPS=3.d-14)
!     Given the lower and upper limits of integration x1 and x2, and given n, this routine returns arrays x(1:n) and w(1:n) of length n, containing the abscissas and weights of the Gauss-Legendre n-point quadrature formula.
INTEGER i,j,m
DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
!     High precision is a good idea for this routine.
 !parece que n debe ser impar
m=(n+1)/2   !The roots are symmetric in the interval, so we only have to ﬁnd half of them.
xm=0.5d0*(x2+x1)
xl=0.5d0*(x2-x1)
do i=1,m  !Loop over the desired roots.
     z=cos(3.141592654d0*(i*1.d0-.25d0)/(n*1.d0+.5d0))
!Starting with the above approximation to the ith root, we enter the main loop of re-
!        ﬁnement by Newton’s method.
1   continue
         p1=1.d0
         p2=0.d0
! Loop up the recurrence relation to get the Leg-endre polynomial evaluated at z.
         do j=1,n
              p3=p2
              p2=p1
              p1=((2.d0*j-1.d0)*z*p2-(j*1.d0-1.d0)*p3)/(j*1.d0)
         enddo
! p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a standard relation involving also p2, the polynomial of one lower order.
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
         z=z1-p1/pp         !Newton’s method.
    if(abs(z-z1).gt.EPS)goto 1
! Scale the root to the desired interval,and put in its symmetric counterpart.
    x(i)=xm-xl*z
    x(n+1-i)=xm+xl*z
    w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)     !Compute the weight and its symmetric counterpart.
    w(n+1-i)=w(i)
enddo
return
END

