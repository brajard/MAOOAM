
! int_comp.f90
!
!> Utility module containing the routines to perform the integration of functions 
!
!> @copyright                                                               
!> 2018 Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!> Most are taken from the Numerical Recipes
!                                                                           
!---------------------------------------------------------------------------!

MODULE int_comp
  USE stoch_params, only: maxint

  PRIVATE

  PUBLIC :: integrate
  
  
CONTAINS
  !> Routine to compute integrals of function from O to #maxint
  !> @param func function to integrate
  !> @param ss result of the integration
  SUBROUTINE integrate(func,ss)
    REAL(KIND=8) :: ss,func,b
    EXTERNAL func
    b=maxint
    ! CALL qromo(func,0.D0,1.D0,ss,midexp)
    CALL qromb(func,0.D0,b,ss)
  END SUBROUTINE integrate

  !> Romberg integration routine
  !> @param func function to integrate
  !> @param a lower limit of the integral
  !> @param b higher limit of the integral
  !> @param func function to integrate
  !> @param ss result of the integration
  SUBROUTINE qromb(func,a,b,ss)
    INTEGER :: JMAX,JMAXP,K,KM
    REAL(KIND=8) :: a,b,func,ss,EPS
    EXTERNAL func
    PARAMETER (EPS=1.D-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
    INTEGER j
    REAL(KIND=8) :: dss,h(JMAXP),s(JMAXP)
    h(1)=1.
    DO j=1,JMAX
       CALL trapzd(func,a,b,s(j),j)
       IF (j.ge.K) THEN
          CALL polint(h(j-KM),s(j-KM),K,0.D0,ss,dss)
          IF (abs(dss).le.EPS*abs(ss)) RETURN
       ENDIF
       s(j+1)=s(j)
       h(j+1)=0.25*h(j)
    ENDDO
    stop 'too many steps in qromb'
  END SUBROUTINE qromb

  !> Romberg integration routine on an open interval
  !> @param a lower limit of the integral
  !> @param b higher limit of the integral
  !> @param func function to integrate
  !> @param ss result of the integration
  !> @param chose routine to perform the integration
  SUBROUTINE qromo(func,a,b,ss,choose)
    INTEGER :: JMAX,JMAXP,K,KM
    REAL(KIND=8) :: a,b,func,ss,EPS
    EXTERNAL func,choose
    PARAMETER (EPS=1.e-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
    INTEGER :: j
    REAL(KIND=8) :: dss,h(JMAXP),s(JMAXP)
    h(1)=1.
    DO j=1,JMAX
       CALL choose(func,a,b,s(j),j)
       IF (j.ge.K) THEN
          call polint(h(j-KM),s(j-KM),K,0.D0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
       ENDIF
       s(j+1)=s(j)
       h(j+1)=h(j)/9.
    ENDDO
    stop 'too many steps in qromo'
  END SUBROUTINE qromo

  !> Polynomial interpolation routine 
  SUBROUTINE polint(xa,ya,n,x,y,dy)
    INTEGER :: n,NMAX
    REAL(KIND=8) :: dy,x,y,xa(n),ya(n)
    PARAMETER (NMAX=10)
    INTEGER :: i,m,ns
    REAL(KIND=8) :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=abs(x-xa(1))
    do i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
    end do
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)stop 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       end do
       if (2*ns.lt.n-m)then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    end do
    return
  END SUBROUTINE polint

  !> Trapezoidal rule integration routine
  !> @param func function to integrate
  !> @param a lower limit of the integral
  !> @param b higher limit of the integral
  !> @param s result of the integration
  !> @param n higher stage of the rule to be computed
  SUBROUTINE trapzd(func,a,b,s,n)
    INTEGER :: n
    REAL(KIND=8) :: a,b,s,func
    EXTERNAL func
    INTEGER :: it,j
    REAL(KIND=8) :: del,sum,tnm,x
    if (n.eq.1) then
       s=0.5*(b-a)*(func(a)+func(b))
    else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do j=1,it
          sum=sum+func(x)
          x=x+del
       end do
       s=0.5*(s+(b-a)*sum/tnm)
    endif
    return
  END SUBROUTINE trapzd

  !> Midpoint rule integration routine
  !> @param func function to integrate
  !> @param a lower limit of the integral
  !> @param b higher limit of the integral
  !> @param s result of the integration
  !> @param n higher stage of the rule to be computed
  SUBROUTINE midpnt(func,a,b,s,n)
    INTEGER :: n
    REAL(KIND=8) :: a,b,s,func
    EXTERNAL func
    INTEGER :: it,j
    REAL(KIND=8) :: ddel,del,sum,tnm,x
    if (n.eq.1) then
       s=(b-a)*func(0.5*(a+b))
    else
       it=3**(n-2)
       tnm=it
       del=(b-a)/(3.*tnm)
       ddel=del+del
       x=a+0.5*del
       sum=0.
       do j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
       end do
       s=(s+(b-a)*sum/tnm)/3.
    endif
    return
  END SUBROUTINE midpnt

  !> Midpoint routine for bb infinite with funk decreasing infinitely rapidly at
  !> infinity
  !> @param funk function to integrate
  !> @param aa lower limit of the integral
  !> @param bb higher limit of the integral
  !> @param s result of the integration
  !> @param n higher stage of the rule to be computed
  SUBROUTINE midexp(funk,aa,bb,s,n)
    INTEGER :: n
    REAL(KIND=8) :: aa,bb,s,funk
    EXTERNAL funk
    INTEGER :: it,j
    REAL(KIND=8) :: ddel,del,sum,tnm,x,func,a,b
    func(x)=funk(-log(x))/x
    b=exp(-aa)
    a=0.
    if (n.eq.1) then
       s=(b-a)*func(0.5*(a+b))
    else
       it=3**(n-2)
       tnm=it
       del=(b-a)/(3.*tnm)
       ddel=del+del
       x=a+0.5*del
       sum=0.
       do j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
       end do
       s=(s+(b-a)*sum/tnm)/3.
    endif
    return
  END SUBROUTINE midexp
END MODULE int_comp
