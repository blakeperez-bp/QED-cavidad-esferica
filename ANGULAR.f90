!Integral de Gaunt: \int d\Omega Y_{l1,m1}^{*}Y_{l2,m2}Y_{l3,m3}

module funciones_angulares
  use pi_module
  implicit none
  private
  public :: I4, cp, dp_minuss, dp_pluss, dpluss, dminuss, I3, II3   ! exportamos solo lo necesario

  contains


    real*8 Function I4(l1,m1,l2,m2,l3,m3,l4,m4)  !!esta no considera que el lpp mpp es conjugado para la integral
      use pi_module
      implicit none
      integer l1, m1, l2, m2, l3, m3, l4, m4
      integer l, Lmmin, Lmmax
      real*8  CG1, CG2, CG3, CG4, fact


      I4 = 0.d0

      Lmmin = max(max(abs(l1-l2),abs(m1+m2)),max(abs(l3-l4),abs(m3+m4)))
      Lmmax = min(l1+l2,l3+l4)

      If(Lmmin.gt.Lmmax) return
      If( (m1+m2+m3+m4) .ne. 0 ) return

      Do l = Lmmin, Lmmax
        fact =                         dsqrt(((2*l1+1)*(2*l2+1))/(4.d0*pi*(2*l+1)))
        fact = fact*((-1.d0)**(m1+m2))*dsqrt(((2*l3+1)*(2*l4+1))/(4.d0*pi*(2*l+1)))
        call clebsch(l1*1.d0,l2*1.d0,l*1.d0,    0.d0,    0.d0,         0.d0,CG1)
        call clebsch(l1*1.d0,l2*1.d0,l*1.d0, m1*1.d0, m2*1.d0, (m1+m2)*1.d0,CG2)
        call clebsch(l3*1.d0,l4*1.d0,l*1.d0,    0.d0,    0.d0,         0.d0,CG3)
        call clebsch(l3*1.d0,l4*1.d0,l*1.d0, m3*1.d0, m4*1.d0,-(m1+m2)*1.d0,CG4)  !tiene la delta adentro
        I4 = I4 + fact*CG1*CG2*CG3*CG4
      End Do

    End


    real*8 Function I3(l1,m1,l2,m2,l3,m3)  !!enteoria esto es para Yl1m1 conjugado
      use pi_module
      implicit none
      integer, intent(in) :: l1,m1,l2,m2,l3,m3
      real*8:: CG1,CG2
      call clebsch(l1*1.d0,l2*1.d0,l3*1.d0,0.d0,0.d0,0.d0,CG1)
      call clebsch(l1*1.d0,l2*1.d0,l3*1.d0,-m1*1.d0,m2*1.d0,-m3*1.d0,CG2)
      I3=((-1.d0)**(m3+m1))*dsqrt(((2*l1+1)*(2*l2+1))/(4.d0*pi*(2*l3+1)))*CG1*CG2
    End
    !!ESTA FUNCION DA BIEN

    real*8 function II3(l1,m1,l2,m2,l3,m3) result(res)  
      use pi_module
      implicit none
      integer, intent(in) :: l1,m1,l2,m2,l3,m3
      res = I4(l1,-m1, 0, 0, l2, m2, l3, m3)*sqrt(4*pi)*(-1)**(m1)  !!ahora aparece que lpp es conjugado
    End



    ! Funcion para evaluar los dlm+-, clm

    real*8 function cp(l,m) result(res)
      integer, intent(in) :: l,m
      res = sqrt( ( (l*1.d0)**2 - (m*1.d0)**2  ) / ( (2*l-1)*(2*l+1)*1.d0  ) )
    End function cp

    real*8 function dp_pluss(l,m) result(res)
      integer, intent(in) :: l,m
      res = sqrt( ( l - m - 1 )*(l - m)*1.d0 / (2*(2*l-1)*(2*l + 1)*1.d0) )
    End function dp_pluss

    real*8 function dpluss(l,m) result(res)
      integer, intent(in) :: l,m
      res = sqrt( ( l + m + 1 )*(l + m + 2)*1.d0 / (2*(2*l+1)*(2*l + 3)*1.d0) )
    end function dpluss


    real*8 function dp_minuss(l,m) result(res)
      integer, intent(in) :: l,m
      res = sqrt( ( l + m - 1 )*(l + m)*1.d0 / (2*(2*l-1)*(2*l+1)*1.d0) )
    End function dp_minuss

    real*8 function dminuss(l,m) result(res)
      integer, intent(in) :: l,m
      res = sqrt( ( l - m + 1 )*(l - m + 2)*1.d0 / (2*(2*l+1)*(2*l + 3)*1.d0) )
    end function dminuss



end module



!!UTIL: comentario sobre coeficientes primados
!! c_{lm} = c'_{l+1,m}
!!dplus_{lm} = dp_plus_{l+1,-m+1}
!!dminus_{lm} = dp_minus_{l+1,-m-1}
!!los no primados estan en funcion de los primados

subroutine evaluateangular1(lpp,mpp,l,m,lp,mp,res1)


  use funciones_angulares
  implicit none
  integer, intent(in) :: lpp, mpp, l, m ,lp, mp  !indices que ingreso
  real*8, intent(out) :: res1 !nro que me va a devolver
  ! real(8) :: cp, dp_minus, dp_plus, I4   !
  ! external :: cp, dp_minus, dp_plus, I4  !
  real*8 ::  c,ccp, dplus,dpplus,  dminus, dpminus, c_p, cp_p, dplus_p,dpplus_p, dminus_p, dpminus_p

  !!como dentro de esto l y m constantes defino:
  c = cp(l+1,m)
  ccp = cp(l,m)
  dplus = dpluss(l,m)
  dpplus = dp_pluss(l,m)
  dminus = dminuss(l,m)
  dpminus = dp_minuss(l,m)

  !!los coeficientes para lp y mp
  c_p = cp(lp+1,mp)
  cp_p = cp(lp,mp)
  dplus_p = dpluss(lp,mp)
  dpplus_p = dp_pluss(lp,mp)
  dminus_p = dminuss(lp,mp)
  dpminus_p = dp_minuss(lp,mp)

  !I4(l1,m1,l2,m2,l3,m3,l4,m4)
  !!1==pp (material izq)
  !!2==1 o 0 (parte que viene del r)
  !!3 == modo fotonico
  !!4 == p (material der)

  res1 = 1.d0 * dminus_p * ( &
        + dplus * l    * I4(lpp, mpp, 1, 0, l+1, m+1, lp+1, mp-1) &
        + dpplus*(l+1) * I4(lpp, mpp, 1, 0, l-1, m+1, lp+1, mp-1) &
        - c     * l    * I4(lpp, mpp, 1, 1, l+1, m  , lp+1, mp-1) &
        + ccp   *(l+1) * I4(lpp, mpp, 1, 1, l-1, m  , lp+1, mp-1) ) &
     + c_p * ( &
        + dminus* l    * I4(lpp, mpp, 1, 1, l+1, m-1, lp+1, mp  ) &
        + dpminus*(l+1)* I4(lpp, mpp, 1, 1, l-1, m-1, lp+1, mp  ) &
        - dplus  * l   * I4(lpp, mpp, 1,-1, l+1, m+1, lp+1, mp  ) &
        - dpplus *(l+1)* I4(lpp, mpp, 1, 1, l-1, m+1, lp+1, mp  ) ) &
     - dplus_p * ( &
        + dminus * l   * I4(lpp, mpp, 1, 0, l+1, m-1, lp+1, mp+1) &
        + dpminus*(l+1)* I4(lpp, mpp, 1, 1, l-1, m-1, lp+1, mp+1) &
        - c      * l   * I4(lpp, mpp, 1,-1, l+1, m  , lp+1, mp+1) &
        + ccp    *(l+1)* I4(lpp, mpp, 1,-1, l-1, m  , lp+1, mp+1) )


End subroutine evaluateangular1


subroutine evaluateangular2(lpp,mpp,l,m,lp,mp, res2)
  use funciones_angulares
  implicit none

  integer, intent(in) :: lpp, mpp, l, m ,lp, mp  !indices que ingreso
  real*8, intent(out) :: res2 !nro que me va a devolver
  ! real(8) :: cp, dp_minus, dp_plus, I4   ! 
  ! external :: cp, dp_minus, dp_plus, I4  ! 
  real*8 ::  c,ccp, dplus,dpplus,  dminus, dpminus, c_p, cp_p, dplus_p,dpplus_p, dminus_p, dpminus_p 

  !!como dentro de esto l y m constantes defino:
  c = cp(l+1,m)
  ccp = cp(l,m)
  dplus = dpluss(l,m)
  dpplus = dp_pluss(l,m)
  dminus = dminuss(l,m)
  dpminus = dp_minuss(l,m)

  !!los coeficientes para lp y mp
  c_p = cp(lp+1,mp)
  cp_p = cp(lp,mp)
  dplus_p = dpluss(lp,mp)
  dpplus_p = dp_pluss(lp,mp)
  dminus_p = dminuss(lp,mp)
  dpminus_p = dp_minuss(lp,mp)


  !I4(l1,m1,l2,m2,l3,m3,l4,m4)
  !!1==pp (material izq)
  !!2==1 o 0 (parte que viene del r)
  !!3 == modo fotonico
  !!4 == p (material der)


  res2 = 1.d0 * dpminus_p * ( &
          - dplus * l    * I4(lpp, mpp, 1, 0, l+1, m+1, lp-1, mp-1) &
          - dpplus*(l+1) * I4(lpp, mpp, 1, 0, l-1, m+1, lp-1, mp-1) &
          + c     * l    * I4(lpp, mpp, 1, 1, l+1, m  , lp-1, mp-1) &
          - ccp   *(l+1) * I4(lpp, mpp, 1, 1, l-1, m  , lp-1, mp-1) ) &
      + cp_p * ( &
          + dminus* l    * I4(lpp, mpp, 1, 1, l+1, m-1, lp-1, mp  ) &
          + dpminus*(l+1)* I4(lpp, mpp, 1, 1, l-1, m-1, lp-1, mp  ) &
          - dplus  * l   * I4(lpp, mpp, 1,-1, l+1, m+1, lp-1, mp  ) &
          - dpplus *(l+1)* I4(lpp, mpp, 1,-1, l-1, m+1, lp+1, mp  ) ) &
      + dpplus_p * ( &
          + dminus * l   * I4(lpp, mpp, 1, 0, l+1, m-1, lp-1, mp+1) &
          + dpminus*(l+1)* I4(lpp, mpp, 1, 0, l-1, m-1, lp-1, mp+1) &
          - c      * l   * I4(lpp, mpp, 1,-1, l+1, m  , lp-1, mp+1) &
          + ccp    *(l+1)* I4(lpp, mpp, 1,-1, l-1, m  , lp-1, mp+1) )


End subroutine evaluateangular2




subroutine factor_angular(lpp,mpp,l,m,lp,mp, res)
  use funciones_angulares
  use pi_module
  implicit none

  integer, intent(in) :: lpp, mpp, l, m ,lp, mp  !indices que ingreso
  real*8, intent(out) :: res!nro que me va a devolver
  real*8 :: sum
  integer :: j, jmax
  

  res = 0.d0

  jmax = int((l + lp + lpp -1)*0.5d0)+1
  do j = 0, jmax
    sum =   c(j)*(mp*sqrt( l*(l+1)*1.d0 - m*(m-1)*1.d0   )* I4(lpp, -mpp, 2*j+1, 1, l, m-1, lp, mp) &
          - m*sqrt( lp*(lp+1)*1.d0 - mp*(mp-1)*1.d0)* I4(lpp, -mpp, 2*j+1, 1, l, m, lp, mp-1) &
          + mp*sqrt(  l*(l+1)*1.d0 - m*(m+1)*1.d0  )* I4(lpp, -mpp, 2*j+1, -1, l, m+1, lp, mp)&
          - m*sqrt( lp*(lp+1)*1.d0 - mp*(mp+1)*1.d0)*I4(lpp, -mpp, 2*j+1, -1, l, m, lp, mp+1))

    res = res + sum
  end do
  res = res*(-1)**mpp  !!para considerar que lo essandwicho con el conjugado

  contains
    real*8 function c(j) result(r)
      integer, intent(in):: j
      r = sqrt(  pi*(4.d0*j + 3.d0) / ( (2.d0*j + 1.d0)*(2.d0*j+2.d0) )   )
    end function c


end subroutine factor_angular






!Integral de Gaunt: \int d\Omega Y_{l1,m1}^{*}Y_{l2,m2}Y_{l3,m3}

! ''''''''''''''''''''''''''''''''''''''''''''''

! !Integrales 2D!Integrales 2D!Integrales 2D!Integrales 2D!Integrales 2D
! !       call AM(Ltot,Mtot,Ltot,Mtot,ele,l1p,l2p,l1,l2,IAM)
! subroutine AM(Lt,Mt,Jt,Nt,l,l1,l2,l3,l4,IAM)
! implicit none
! integer Lt,Mt,Jt,Nt,l,l1,l2,l3,l4
! integer m,m1,m3,minm1,maxm1,minm3,maxm3
! real*8 CG1,CG2,CG3,CG4,CG5,CG6,IAM
! real*8 Suma1,Suma2,Suma3
! IAM=0.d0
! call clebsch(l3*1.d0,l*1.d0,l1*1.d0,0.d0,0.d0,0.d0,CG1)
! call clebsch(l4*1.d0,l*1.d0,l2*1.d0,0.d0,0.d0,0.d0,CG2)
! IAM=CG1*CG2*sqrt(((2.d0*l3+1.d0)*(2.d0*l4+1.d0))/((2.d0*l1+1.d0)*(2.d0*l2+1.d0)))
! If((IAM.ne.0.d0).and.((Lt.eq.Jt).and.(Mt.eq.Nt)))Then
! minm1=max(-l1,Mt-l2)
! maxm1=min(l1,Mt+l2)
! minm3=max(-l3,Mt-l4)
! maxm3=min(l3,Mt+l4)
! Suma1=0.d0
!       DO m3=minm3,maxm3
!       Suma2=0.d0
!             DO m1=minm1,maxm1
!             Suma3=0.d0
!                   Do m=-l,l
!                   call clebsch(l3*1.d0,l*1.d0,l1*1.d0,-m3*1.d0,m*1.d0,-m1*1.d0,CG3)
!                   call clebsch(l4*1.d0,l*1.d0,l2*1.d0,(m3-Mt)*1.d0,-m*1.d0,(m1-Mt)*1.d0,CG4)
!                   Suma3=Suma3+((-1.d0)**m)*CG3*CG4
!                   End Do
!             call clebsch(l1*1.d0,l2*1.d0,Lt*1.d0,m1*1.d0,(Mt-m1)*1.d0,Mt*1.d0,CG5)
!             Suma2=Suma2+Suma3*CG5
!             End Do
!       call clebsch(l3*1.d0,l4*1.d0,Lt*1.d0,m3*1.d0,(Mt-m3)*1.d0,Mt*1.d0,CG6)
!       Suma1=Suma1+Suma2*CG6
!       END DO
! IAM=IAM*Suma1
! End If
! end subroutine AM

! subroutine AMolecular(Lt,Mt,Jt,Nt,l,l1,l2,l1p,l2p,IAM)
! implicit none
! integer Lt,Mt,Jt,Nt,l,l1,l2,l1p,l2p
! integer m,m1,m1p,minm1,maxm1,minm3,maxm3
! real*8 CG1,CG2,CG3,CG4,CG5,CG6,IAM
! real*8 Suma1,Suma2,Suma3
! real*8, external :: I3
! IAM=0.d0
! If((Lt.eq.Jt).and.(Mt.eq.Nt))Then
! minm1=max(-l1,Mt-l2)
! maxm1=min(l1,Mt+l2)
! minm3=max(-l1p,Mt-l2p)
! maxm3=min(l1p,Mt+l2p)
! Suma1=0.d0
!       DO m1p=minm3,maxm3
!       Suma2=0.d0
!             DO m1=minm1,maxm1
!             Suma3=0.d0
!                   Do m=-l,l
! 		   CG3=I3(l,m,l1p,-m1p,l1,m1)
! 		   CG4=I3(l,-m,l2p,m1p-Mt,l2,Mt-m1)
!                   Suma3=Suma3+((-1.d0)**m)*CG3*CG4
!                   End Do
!             call clebsch(l1*1.d0,l2*1.d0,Lt*1.d0,m1*1.d0,(Mt-m1)*1.d0,Mt*1.d0,CG5)
!             Suma2=Suma2+Suma3*CG5
!             End Do
!       call clebsch(l1p*1.d0,l2p*1.d0,Lt*1.d0,m1p*1.d0,(Mt-m1p)*1.d0,Mt*1.d0,CG6)
!       Suma1=Suma1+Suma2*CG6
!       END DO
! IAM=Suma1*2.d0*3.14159265359d0*((-1.d0)**Mt)
! End If
! end subroutine AMolecular


! !integrales 1D!integrales 1D!integrales 1D!integrales 1D!integrales 1D

! subroutine ProdescD(Lt,Mt,l1p,l2p,l1,l2,PD) !!FUNCIONA
! implicit none
! integer Mt,Lt,l1p,l2p,l1,l2
! integer m,m1,m1p,minl1l2,maxl1l2
! real*8 CG1,PD,Suma
! If((l1p.eq.l1).and.(l2p.eq.l2))Then
! PD=1.d0
! Else
! PD=0.d0
! End If
! !   If((l1p.eq.l1).and.(l2p.eq.l2))Then
! !   minl1l2=max(-l1,Mt-l2)
! !   maxl1l2=min(l1,Mt+l2)
! !   Suma=0.d0
! !         DO m1=minl1l2,maxl1l2
! !         call clebsch(l1*1.d0,l2*1.d0,Lt*1.d0,m1*1.d0,(Mt-m1)*1.d0,Mt*1.d0,CG1)
! !         Suma=Suma+CG1*CG1
! !         End Do
! !   PD=Suma
! !   Else
! !   PD=0.d0
! !   End If
! end subroutine ProdescD

! subroutine ProdescI(Lt,Mt,l1p,l2p,l1,l2,PI) !!FUNCIONA
! implicit none
! integer Mt,Lt,l1p,l2p,l1,l2
! integer m,m1,m1p,minl1l2,maxl1l2,minl1pl2p,maxl1pl2p
! real*8 CG1,CG2,PI,Suma1,Suma2
!       If((l2p.eq.l1).and.(l1p.eq.l2))Then
!       minl1pl2p=max(-l2,Mt-l1)
!       maxl1pl2p=min(l2,Mt+l1)
!       minl1l2=max(-l1,Mt-l2)
!       maxl1l2=min(l1,Mt+l2)
!       Suma1=0.d0
!             DO m1p=minl1pl2p,maxl1pl2p
!             Suma2=0.d0
!                   DO m1=minl1l2,maxl1l2
!                   CG1=0.d0
!                         If(m1p.eq.(Mt-m1))Then
!                         call clebsch(l1*1.d0,l2*1.d0,Lt*1.d0,m1*1.d0,(Mt-m1)*1.d0,Mt*1.d0,CG1)
!                         End If
!                   Suma2=Suma2+CG1
!                   End Do
!             call clebsch(l2*1.d0,l1*1.d0,Lt*1.d0,m1p*1.d0,(Mt-m1p)*1.d0,Mt*1.d0,CG2)
!             Suma1=Suma1+CG2*Suma2
!             End Do
!       PI=Suma1
!       Else
!       PI=0.d0
!       End If
! end subroutine ProdescI

! !integrales 1D!integrales 1D!integrales 1D!integrales 1D!integrales 1D



! !MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ!MATRIZ

! !ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!
! !ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!
! !ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!
! !ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!


! subroutine Pnortogonal(Lp,Mp,lap,lbp,L,M,la,lb,IA_3)
! use pi_module
! use complex_module
! implicit none
! integer Lp,Mp,lap,lbp,L,M,la,lb,l1,l2
! integer ma,mb,map,mbp,m1,m2
! real*8 CG1,CG2,Pno,Suma1,Suma2
! real*8, external :: I3,A_exp,A_sinm1,A_ctag
! complex*16, external :: a1_l_m,a2_l_m,a3_l_m,a3a_l_m,a3b_l_m,a1a_l_m,a2a_l_m,a1b_l_m,a2b_l_m,a1c_l_m,a2c_l_m
! complex*16 IA1,IA2,IA3,IA4,IA5,IA6,IA7,IA8,IA9,IA10,IA11,IA12
! complex*16 IA_3(9),I3A,I3Ap1,I3B,I3Bp1
! complex*16 prueba1,prueba2,prueba3,prueba4,prueba5,prueba6		!BORRAR
! complex*16 prueba7,prueba8,prueba9,prueba10,prueba11,prueba12		!BORRAR
! complex*16 a1_c1,a2_c1,a3_c1,a3a_c1,a3b_c1,a1a_c1
! complex*16 a1_c2,a2_c2,a3_c2,a3a_c2,a3b_c2,a1a_c2
! complex*16 a2a_c1,a1b_c1,a2b_c1,a1c_c1,a2c_c1
! complex*16 a2a_c2,a1b_c2,a2b_c2,a1c_c2,a2c_c2
! IA_3=c0

! ! prueba1=c0
! ! prueba2=c0
! ! prueba3=c0
! ! prueba4=c0
! ! prueba5=c0
! ! prueba6=c0
! ! prueba7=c0
! ! prueba8=c0
! ! prueba9=c0
! ! prueba10=c0
! ! prueba11=c0
! ! prueba12=c0

! ! write(*,*) lap,lbp,la,lb

! Do map=max(-lap,Mp-lbp),Min(lap,Mp+lbp)
! mbp=Mp-map
! call clebsch(lap*1.d0,lbp*1.d0,Lp*1.d0,map*1.d0,mbp*1.d0,Mp*1.d0,CG1)
!   Do ma=max(-la,M-lb),Min(la,M+lb)
!   mb=M-ma
!   call clebsch(la*1.d0,lb*1.d0,L*1.d0,ma*1.d0,mb*1.d0,M*1.d0,CG2)


!     Do l1=int(abs(la-lap)*1.001d0),(la+lap)
!       Do m1=map-ma-1,map-ma			!-l1,l1
! !       Do m1=-l1,l1


! 	Do l2=int(abs(lb-lbp)*1.001d0),(lb+lbp)
! 	  Do m2=ma-map-1,ma-map				!-l2,l2
! ! 	  Do m2=-l2,l2

! 	  If(int(abs(m1)*1.001d0).gt.l1)Then
! 	  goto	20
! 	  End If
! 	  If(int(abs(m2)*1.001d0).gt.l2)Then
! 	  goto	20
! 	  End If
! 	  If(CG1*CG2.eq.0.d0)Then
! 	  goto	20
! 	  End If
! ! If(int(abs(map-ma)*1.001d0).le.l1)Then
! ! m1=map-ma
! ! I3A=c1*I3(lap,map,la,ma,l1,m1)
! ! I3Ap1=c0
! ! ElseIf(int(abs(map-ma-1)*1.001d0).le.l1)Then
! ! I3Ap1=c1*I3(lap,map,la,ma+1,l1,m1)
! ! I3A=c0
! ! End If

! I3A=c1*I3(lap,map,la,ma,l1,m1)
! I3Ap1=c1*I3(lap,map,la,ma+1,l1,m1)*sqrt(1.d0*(la-ma)*(1+la+ma))
! I3B=c1*I3(lbp,Mp-map,lb,Mp-ma,l2,m2)
! I3Bp1=c1*I3(lbp,Mp-map,lb,Mp-ma+1,l2,m2)*sqrt(1.d0*(lb-M+ma)*(1+lb+M-ma))
! ! 
! ! write(*,*) "---------------"
! ! write(*,*) map,ma
! ! write(*,*) l1,m1
! ! write(*,*) l2,m2
! ! write(*,*) "---------------"
! a1_c1=a1_l_m(l1,m1)
! a1_c2=a1_l_m(l2,m2)

! a2_c1=a2_l_m(l1,m1)
! a2_c2=a2_l_m(l2,m2)

! a3_c1=a3_l_m(l1,m1)
! a3_c2=a3_l_m(l2,m2)

! a3a_c1=a3a_l_m(l1,m1)
! a3a_c2=a3a_l_m(l2,m2)

! a3b_c1=a3b_l_m(l1,m1)
! a3b_c2=a3b_l_m(l2,m2)

! a1a_c1=a1a_l_m(l1,m1)
! a1a_c2=a1a_l_m(l2,m2)

! a2a_c1=a2a_l_m(l1,m1)
! a2a_c2=a2a_l_m(l2,m2)


! a2b_c1=a2b_l_m(l1,m1)
! a2b_c2=a2b_l_m(l2,m2)

! a1b_c1=a1b_l_m(l1,m1)
! a1b_c2=a1b_l_m(l2,m2)

! a1c_c1=a1c_l_m(l1,m1)
! a1c_c2=a1c_l_m(l2,m2)

! a2c_c1=a2c_l_m(l1,m1)
! a2c_c2=a2c_l_m(l2,m2)


! !111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
! IA1=a1_c1*a1_c2*I3A*I3B
! IA2=a2_c1*a2_c2*I3A*I3B
! IA3=a3_c1*a3_c2*I3A*I3B
! ! prueba1=prueba1+CG1*CG2*IA1
! ! prueba2=prueba2+CG1*CG2*IA2
! ! prueba3=prueba3+CG1*CG2*IA3
! IA_3(1)=IA_3(1)+CG1*CG2*(IA1+IA2+IA3)
! !111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111

! !222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
! IA1= a1_c1*a1a_c2*(M-ma)*I3A*I3B
! IA2= a1_c1*a1b_c2*I3A*I3Bp1
! IA3= a2_c1*a2a_c2*(M-ma)*I3A*I3B
! IA4= a2_c1*a2b_c2*I3A*I3Bp1
! IA5=-a3_c1*a3a_c2*(M-ma)*I3A*I3B
! IA6=-a3_c1*a3b_c2*I3A*I3Bp1
! ! prueba1=prueba1+CG1*CG2*IA1
! ! prueba2=prueba2+CG1*CG2*IA2
! ! prueba3=prueba3+CG1*CG2*IA3
! ! prueba4=prueba4+CG1*CG2*IA4
! ! prueba5=prueba5+CG1*CG2*IA5
! ! prueba6=prueba6+CG1*CG2*IA6
! IA_3(2)=IA_3(2)+CG1*CG2*(IA1+IA2+IA3+IA4+IA5+IA6)
! !222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222

! !333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
! IA1=-ci*(M-ma)*a1_c1*a1c_c2*I3A*I3B
! IA2= ci*(M-ma)*a2_c1*a2c_c2*I3A*I3B
! ! prueba1=prueba1+CG1*CG2*IA1
! ! prueba2=prueba2+CG1*CG2*IA2
! IA_3(3)=IA_3(3)+CG1*CG2*(IA1+IA2)
! !333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333

! !444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
! IA1= a1a_c1*a1_c2*ma*I3A*I3B
! IA2= a1b_c1*a1_c2*I3Ap1*I3B
! IA3= a2a_c1*a2_c2*ma*I3A*I3B
! IA4= a2b_c1*a2_c2*I3Ap1*I3B
! IA5=-a3a_c1*a3_c2*ma*I3A*I3B
! IA6=-a3b_c1*a3_c2*I3Ap1*I3B
! ! prueba1=prueba1+CG1*CG2*IA1
! ! prueba2=prueba2+CG1*CG2*IA2
! ! prueba3=prueba3+CG1*CG2*IA3
! ! prueba4=prueba4+CG1*CG2*IA4
! ! prueba5=prueba5+CG1*CG2*IA5
! ! prueba6=prueba6+CG1*CG2*IA6
! IA_3(4)=IA_3(4)+CG1*CG2*(IA1+IA2+IA3+IA4+IA5+IA6)
! !444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444

! !555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
! IA1=a1a_c1*a1a_c2*ma*(M-ma)*I3A*I3B
! IA2=a1b_c1*a1a_c2*(M-ma)*I3Ap1*I3B
! IA3=a1a_c1*a1b_c2*ma*I3A*I3Bp1
! IA4=a1b_c1*a1b_c2*I3Ap1*I3Bp1

! IA5=a2a_c1*a2a_c2*ma*(M-ma)*I3A*I3B
! IA6=a2b_c1*a2a_c2*(M-ma)*I3Ap1*I3B
! IA7=a2a_c1*a2b_c2*ma*I3A*I3Bp1
! IA8=a2b_c1*a2b_c2*I3Ap1*I3Bp1

! IA9=a3a_c1*a3a_c2*ma*(M-ma)*I3A*I3B
! IA10=a3b_c1*a3a_c2*(M-ma)*I3Ap1*I3B
! IA11=a3a_c1*a3b_c2*ma*I3A*I3Bp1
! IA12=a3b_c1*a3b_c2*I3Ap1*I3Bp1

! ! prueba1=prueba1+CG1*CG2*IA1
! ! prueba2=prueba2+CG1*CG2*IA2
! ! prueba3=prueba3+CG1*CG2*IA3
! ! prueba4=prueba4+CG1*CG2*IA4
! ! prueba5=prueba5+CG1*CG2*IA5
! ! prueba6=prueba6+CG1*CG2*IA6
! ! prueba7=prueba7+CG1*CG2*IA7
! ! prueba8=prueba8+CG1*CG2*IA8
! ! prueba9=prueba9+CG1*CG2*IA9
! ! prueba10=prueba10+CG1*CG2*IA10
! ! prueba11=prueba11+CG1*CG2*IA11
! ! prueba12=prueba12+CG1*CG2*IA12

! IA_3(5)=IA_3(5)+CG1*CG2*(IA1+IA2+IA3+IA4+IA5+IA6+IA7+IA8+IA9+IA10+IA11+IA12)
! !555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555

! !666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666
! IA1= ci*a2a_c1*a2c_c2*ma*(M-ma)*I3A*I3B
! IA2= ci*a2b_c1*a2c_c2*(M-ma)*I3Ap1*I3B
! IA3=-ci*a1a_c1*a1c_c2*ma*(M-ma)*I3A*I3B
! IA4=-ci*a1b_c1*a1c_c2*(M-ma)*I3Ap1*I3B

! ! prueba1=prueba1+CG1*CG2*IA1
! ! prueba2=prueba2+CG1*CG2*IA2
! ! prueba3=prueba3+CG1*CG2*IA3
! ! prueba4=prueba4+CG1*CG2*IA4

! IA_3(6)=IA_3(6)+CG1*CG2*(IA1+IA2+IA3+IA4)
! !666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666

! !777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777
! IA1=-ci*ma*a1c_c1*a1_c2*I3A*I3B
! IA2= ci*ma*a2c_c1*a2_c2*I3A*I3B

! ! prueba1=prueba1+CG1*CG2*IA1
! ! prueba2=prueba2+CG1*CG2*IA2

! IA_3(7)=IA_3(7)+CG1*CG2*(IA1+IA2)
! !777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777

! !888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
! IA1= ci*a2c_c1*a2a_c2*ma*(M-ma)*I3A*I3B
! IA2= ci*a2c_c1*a2b_c2*ma*I3A*I3Bp1
! IA3=-ci*a1c_c1*a1a_c2*ma*(M-ma)*I3A*I3B
! IA4=-ci*a1c_c1*a1b_c2*ma*I3A*I3Bp1

! ! prueba1=prueba1+CG1*CG2*IA1
! ! prueba2=prueba2+CG1*CG2*IA2
! ! prueba3=prueba3+CG1*CG2*IA3
! ! prueba4=prueba4+CG1*CG2*IA4

! IA_3(8)=IA_3(8)+CG1*CG2*(IA1+IA2+IA3+IA4)
! !888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

! !999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999
! IA1=-ma*(M-ma)*a1c_c1*a1c_c2*I3A*I3B
! IA2=-ma*(M-ma)*a2c_c1*a2c_c2*I3A*I3B

! ! prueba1=prueba1+CG1*CG2*IA1
! ! prueba2=prueba2+CG1*CG2*IA2

! IA_3(9)=IA_3(9)+CG1*CG2*(IA1+IA2)
! !999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999

! 20 continue

! 	  End Do
! 	End Do
!       End Do
!     End Do
!   End Do
! End Do

! ! write(*,*) prueba1!+prueba2!+prueba3+prueba4
! ! write(*,*) prueba2
! ! write(*,*) prueba3+prueba4
! ! write(*,*) prueba4
! ! write(*,*) prueba5+prueba6+prueba7+prueba8
! ! write(*,*) prueba6
! ! write(*,*) prueba7
! ! write(*,*) prueba8
! ! write(*,*) prueba9+prueba10+prueba11+prueba12
! ! write(*,*) prueba10
! ! write(*,*) prueba11
! ! write(*,*) prueba12
! ! write(*,*) prueba1+prueba2!+prueba3+prueba4
! ! +prueba5+prueba6+prueba7+prueba8+prueba9+prueba10+prueba11+prueba12

! end subroutine Pnortogonal


! subroutine Gradiente_xyz(Lp,Mp,lap,lbp,L,M,la,lb,g_xyz)
! use pi_module
! use complex_module
! implicit none
! integer Lp,Mp,lap,lbp,L,M,la,lb,le
! integer ma,mb,map,mbp,me
! real*8 CG1,CG2,Pno,Suma1,Suma2
! real*8, external :: I3
! complex*16, external :: a1_l_m,a2_l_m,a3_l_m,a3a_l_m,a3b_l_m,a1a_l_m,a2a_l_m,a1b_l_m,a2b_l_m,a1c_l_m,a2c_l_m
! complex*16 g_xA,g_xB,g_xC,g_yD,g_yE,g_yF,g_zG,g_zH
! complex*16 g_xyz(8)
! g_xyz=c0
! g_xA=c0
! g_xB=c0
! g_xC=c0
! g_yD=c0
! g_yE=c0
! g_yF=c0
! g_zG=c0
! g_zH=c0

! If(lb.ne.lbp)Then
! goto 21
! End If

! Do map=max(-lap,Mp-lbp),Min(lap,Mp+lbp)
! mbp=Mp-map
! call clebsch(lap*1.d0,lbp*1.d0,Lp*1.d0,map*1.d0,mbp*1.d0,Mp*1.d0,CG1)
! ! write(*,*) "limite map",max(-lap,Mp-lbp),Min(lap,Mp+lbp)

!   Do ma=max(-la,M-lb),Min(la,M+lb)
! ! write(*,*) "limite ma",max(-la,M-lb),Min(la,M+lb)
!     If(ma.ne.(M-Mp+map))Then
! !     write(*,*) "Actúa delta(Mp-map,M-ma)"
!     goto 19
!     End If

!   mb=M-ma
!   call clebsch(la*1.d0,lb*1.d0,L*1.d0,ma*1.d0,mb*1.d0,M*1.d0,CG2)
! ! write(*,*) "ma",ma,"map",map
! ! write(*,*) "clebsch",CG1,CG2
! ! pause
!     Do le=int(abs(la-lap)*1.001d0),(la+lap)
!       Do me=-le,le         !map-ma-1,map-ma			!



! 	  If(int(abs(me)*1.001d0).gt.le)Then
! 	  goto	20
! 	  End If

! 	  If(CG1*CG2.eq.0.d0)Then
! 	  goto	20
! 	  End If

! g_xA=g_xA+CG1*CG2*c1*a1_l_m(le,me)*I3(lap,map,la,ma,le,me)

! g_xB=g_xB+CG1*CG2*c1*a1a_l_m(le,me)*ma*I3(lap,map,la,ma,le,me)
! g_xB=g_xB+CG1*CG2*c1*a1b_l_m(le,me)*sqrt(1.d0*(la-ma)*(1+la+ma))*I3(lap,map,la,ma+1,le,me)

! g_xC=g_xC+CG1*CG2*c1*a1c_l_m(le,me)*ci*ma*I3(lap,map,la,ma,le,me)

! g_yD=g_yD+CG1*CG2*c1*a2_l_m(le,me)*I3(lap,map,la,ma,le,me)

! g_yE=g_yE+CG1*CG2*c1*a2a_l_m(le,me)*ma*I3(lap,map,la,ma,le,me)
! g_yE=g_yE+CG1*CG2*c1*a2b_l_m(le,me)*sqrt(1.d0*(la-ma)*(1+la+ma))*I3(lap,map,la,ma+1,le,me)

! g_yF=g_yF+CG1*CG2*c1*a2c_l_m(le,me)*ci*ma*I3(lap,map,la,ma,le,me)

! g_zG=g_zG+CG1*CG2*c1*a3_l_m(le,me)*I3(lap,map,la,ma,le,me)

! ! write(*,*) "le",le,"me",me
! ! write(*,*) a3_l_m(le,me),I3(lap,map,la,ma,le,me)
! ! write(*,*) g_zG,CG1*CG2*c1*a3_l_m(le,me)*I3(lap,map,la,ma,le,me)

! g_zH=g_zH+CG1*CG2*c1*a3_l_m(le,me)*ma*I3(lap,map,la,ma,le,me)
! g_zH=g_zH+CG1*CG2*c1*a3b_l_m(le,me)*sqrt(1.d0*(la-ma)*(1+la+ma))*I3(lap,map,la,ma+1,le,me)


! 20 continue


!       End Do
!     End Do

! 19 continue

!   End Do
! End Do

! g_xyz(1)=g_xA
! g_xyz(2)=g_xB
! g_xyz(3)=g_xC
! g_xyz(4)=g_yD
! g_xyz(5)=g_yE
! g_xyz(6)=g_yF
! g_xyz(7)=g_zG
! g_xyz(8)=g_zH

! 21 continue
! end subroutine Gradiente_xyz


! subroutine Gradiente_x1(Lp,Mp,lap,lbp,L,M,la,lb,g_x)
! use pi_module
! use complex_module
! implicit none
! integer Lp,Mp,lap,lbp,L,M,la,lb,le
! integer ma,mb,map,mbp,me
! real*8 CG1,CG2,Pno,Suma1,Suma2
! real*8, external :: I3
! complex*16, external :: a1_l_m,a2_l_m,a3_l_m,a3a_l_m,a3b_l_m,a1a_l_m,a2a_l_m,a1b_l_m,a2b_l_m,a1c_l_m,a2c_l_m
! complex*16 g_A,g_B,g_C
! complex*16 g_x(3)
! g_x=c0
! g_A=c0
! g_B=c0
! g_C=c0

! If(lb.ne.lbp)Then
! goto 21
! End If

! Do map=max(-lap,Mp-lbp),Min(lap,Mp+lbp)
! mbp=Mp-map
! call clebsch(lap*1.d0,lbp*1.d0,Lp*1.d0,map*1.d0,mbp*1.d0,Mp*1.d0,CG1)

!   Do ma=max(-la,M-lb),Min(la,M+lb)
!     If(ma.ne.(M-Mp+map))Then
!     goto 19
!     End If

!   mb=M-ma
!   call clebsch(la*1.d0,lb*1.d0,L*1.d0,ma*1.d0,mb*1.d0,M*1.d0,CG2)

!     Do le=int(abs(la-lap)*1.001d0),(la+lap)
!       Do me=-le,le         !map-ma-1,map-ma			!



! 	  If(int(abs(me)*1.001d0).gt.le)Then
! 	  goto	20
! 	  End If

! 	  If(CG1*CG2.eq.0.d0)Then
! 	  goto	20
! 	  End If

! g_A=g_A+CG1*CG2*c1*a1_l_m(le,me)*I3(lap,map,la,ma,le,me)

! g_B=g_B+CG1*CG2*c1*a1a_l_m(le,me)*ma*I3(lap,map,la,ma,le,me)
! g_B=g_B+CG1*CG2*c1*a1b_l_m(le,me)*sqrt(1.d0*(la-ma)*(1+la+ma))*I3(lap,map,la,ma+1,le,me)

! g_C=g_C+CG1*CG2*c1*a1c_l_m(le,me)*ci*ma*I3(lap,map,la,ma,le,me)

! 20 continue


!       End Do
!     End Do

! 19 continue

!   End Do
! End Do

! g_x(1)=g_A
! g_x(2)=g_B
! g_x(3)=g_C

! 21 continue
! end subroutine Gradiente_x1

! subroutine Gradiente_y1(Lp,Mp,lap,lbp,L,M,la,lb,g_y)
! use pi_module
! use complex_module
! implicit none
! integer Lp,Mp,lap,lbp,L,M,la,lb,le
! integer ma,mb,map,mbp,me
! real*8 CG1,CG2,Pno,Suma1,Suma2
! real*8, external :: I3
! complex*16, external :: a1_l_m,a2_l_m,a3_l_m,a3a_l_m,a3b_l_m,a1a_l_m,a2a_l_m,a1b_l_m,a2b_l_m,a1c_l_m,a2c_l_m
! complex*16 g_D,g_E,g_F
! complex*16 g_y(3)
! g_y=c0
! g_D=c0
! g_E=c0
! g_F=c0

! If(lb.ne.lbp)Then
! goto 21
! End If

! Do map=max(-lap,Mp-lbp),Min(lap,Mp+lbp)
! mbp=Mp-map
! call clebsch(lap*1.d0,lbp*1.d0,Lp*1.d0,map*1.d0,mbp*1.d0,Mp*1.d0,CG1)

!   Do ma=max(-la,M-lb),Min(la,M+lb)
!     If(ma.ne.(M-Mp+map))Then
!     goto 19
!     End If

!   mb=M-ma
!   call clebsch(la*1.d0,lb*1.d0,L*1.d0,ma*1.d0,mb*1.d0,M*1.d0,CG2)

!     Do le=int(abs(la-lap)*1.001d0),(la+lap)
!       Do me=-le,le         !map-ma-1,map-ma			!



! 	  If(int(abs(me)*1.001d0).gt.le)Then
! 	  goto	20
! 	  End If

! 	  If(CG1*CG2.eq.0.d0)Then
! 	  goto	20
! 	  End If

! g_D=g_D+CG1*CG2*c1*a2_l_m(le,me)*I3(lap,map,la,ma,le,me)

! g_E=g_E+CG1*CG2*c1*a2a_l_m(le,me)*ma*I3(lap,map,la,ma,le,me)
! g_E=g_E+CG1*CG2*c1*a2b_l_m(le,me)*sqrt(1.d0*(la-ma)*(1+la+ma))*I3(lap,map,la,ma+1,le,me)

! g_F=g_F+CG1*CG2*c1*a2c_l_m(le,me)*ci*ma*I3(lap,map,la,ma,le,me)



! 20 continue


!       End Do
!     End Do

! 19 continue

!   End Do
! End Do

! g_y(1)=g_D
! g_y(2)=g_E
! g_y(3)=g_F
! 21 continue
! end subroutine Gradiente_y1

! subroutine Gradiente_z1(Lp,Mp,lap,lbp,L,M,la,lb,g_z)
! use pi_module
! use complex_module
! implicit none
! integer Lp,Mp,lap,lbp,L,M,la,lb,le
! integer ma,mb,map,mbp,me
! real*8 CG1,CG2,Pno,Suma1,Suma2
! real*8, external :: I3
! complex*16, external :: a1_l_m,a2_l_m,a3_l_m,a3a_l_m,a3b_l_m,a1a_l_m,a2a_l_m,a1b_l_m,a2b_l_m,a1c_l_m,a2c_l_m
! complex*16 g_G,g_H
! complex*16 g_z(2)
! g_z=c0


! ! write(*,*) "Gradiente_z1"

! If(lb.ne.lbp)Then
! goto 21
! End If

! Do map=max(-lap,Mp-lbp),Min(lap,Mp+lbp)
!   mbp=Mp-map
!   call clebsch(lap*1.d0,lbp*1.d0,Lp*1.d0,map*1.d0,mbp*1.d0,Mp*1.d0,CG1)

!   Do ma=max(-la,M-lb),Min(la,M+lb)
!     mb=M-ma
!     call clebsch(la*1.d0,lb*1.d0,L*1.d0,ma*1.d0,mb*1.d0,M*1.d0,CG2)

! !       write(*,*) map,ma,mbp,mb,CG1*CG2
    
!     g_G=c0
! g_H=c0
!     If(mb.ne.mbp)then
!       goto 19
!     End IF
    

! ! write(*,*) lap,lbp,la,lb,map,ma,mbp,mb,CG1*CG2


! !     If(CG1*CG2.eq.0.d0)Then
! !     goto	19
! !     End If

!     Do le=int(abs(la-lap)*1.001d0),(la+lap)
!       Do me=-le,le         !map-ma-1,map-ma			!

! ! 	If(int(abs(me)*1.001d0).gt.le)Then
! ! 	  goto	20
! ! 	End If

! 	g_H=g_H+(c1*a3_l_m(le,me)*I3(lap,map,la,ma,le,me))

! 	g_G=g_G+(c1*a3_l_m(le,me)*ma*I3(lap,map,la,ma,le,me))
! 	g_G=g_G+(c1*a3b_l_m(le,me)*sqrt(1.d0*(la-ma)*(1+la+ma))*I3(lap,map,la,ma+1,le,me))

! ! 	20 continue

!       End Do
!     End Do
! !         write(*,*) g_G!CG1*CG2
! !     write(*,*) g_H!CG1*CG2
! g_z(1)=g_z(1)+(g_G*CG1*CG2)
!  g_z(2)=g_z(2)+(g_H*CG1*CG2)
! !  

    
!     19 continue

!   End Do
! End Do

! ! write(*,*) "sale",g_z

! 21 continue
! end subroutine Gradiente_z1

! subroutine GradienteEd_z1(Lp,Mp,lap,lbp,L,M,la,lb,g_z)
! ! empleando las fórmulas de la página 85 del libro de Edmonds Corregido por que la deficinión de los esféricos es diferente
! use pi_module
! use complex_module
! implicit none
! integer Lp,Mp,lap,lbp,L,M,la,lb,le
! integer ma,mb,map,mbp,me
! real*8 CG1,CG2,Pno,Suma1,Suma2
! real*8, external :: I3
! complex*16, external :: a1_l_m,a2_l_m,a3_l_m,a3a_l_m,a3b_l_m,a1a_l_m,a2a_l_m,a1b_l_m,a2b_l_m,a1c_l_m,a2c_l_m
! complex*16 g_G,g_H
! complex*16 g_z(2)
! g_z=c0
! g_G=c0
! g_H=c0
! !     write(*,*)"Edmonds",lap,lbp,la,lb

! If(lb.ne.lbp)Then
! goto 21
! End If

! Do map=max(-lap,Mp-lbp),Min(lap,Mp+lbp)
!   mbp=Mp-map
!   call clebsch(lap*1.d0,lbp*1.d0,Lp*1.d0,map*1.d0,mbp*1.d0,Mp*1.d0,CG1)

!   Do ma=max(-la,M-lb),Min(la,M+lb)
!     mb=M-ma
!     call clebsch(la*1.d0,lb*1.d0,L*1.d0,ma*1.d0,mb*1.d0,M*1.d0,CG2)

!     If(mb.ne.mbp)then
!     goto 19
!     End IF
    
!     g_G=c0
! !seno*partial( Ylm  )
!     Suma1=1.d0*((2*la+1)*(2*la+3))
!     Suma1=1.d0*(((la-ma+1)*(la+ma+1))/Suma1)
!     Suma1=la*(Suma1**0.5d0)
!     If((lap.eq.(la+1)).and.(map.eq.ma))Then
!       g_G=g_G+c1*Suma1
!     End IF

!     Suma1=1.d0*((2*la+1)*(2*la-1))
!     Suma1=1.d0*(((la-ma)*(la+ma))/Suma1)
!     Suma1=(la+1)*(Suma1**0.5d0)
!     If((lap.eq.(la-1)).and.(map.eq.ma))Then
!       g_G=g_G-c1*Suma1
!     End IF
    
!     g_z(1)=g_z(1)+(g_G*CG1*CG2)

    
! !coseno* Ylm
!     g_H=c0

!     Suma1=1.d0*((2*la+1)*(2*la+3))
!     Suma1=1.d0*(((la-ma+1)*(la+ma+1))/Suma1)
!     Suma1=(Suma1**0.5d0)
!     If((lap.eq.(la+1)).and.(map.eq.ma))Then
!       g_H=g_H+c1*Suma1
!     End IF
    
!     Suma1=1.d0*((2*la+1)*(2*la-1))
!     Suma1=1.d0*(((la-ma)*(la+ma))/Suma1)
!     Suma1=(Suma1**0.5d0)
!     If((lap.eq.(la-1)).and.(map.eq.ma))Then
!       g_H=g_H+c1*Suma1
!     End IF

!     g_z(2)=g_z(2)+(g_H*CG1*CG2)
    
!     19  continue
!   End Do
! End Do

! 21 continue
! end subroutine GradienteEd_z1

! subroutine GradienteEd_z2(Lp,Mp,lap,lbp,L,M,la,lb,g_z)
! ! empleando las fórmulas de la página 85 del libro de Edmonds Corregido por que la deficinión de los esféricos es diferente
! use pi_module
! use complex_module
! implicit none
! integer Lp,Mp,lap,lbp,L,M,la,lb,le
! integer ma,mb,map,mbp,me
! real*8 CG1,CG2,Pno,Suma1,Suma2
! real*8, external :: I3
! complex*16, external :: a1_l_m,a2_l_m,a3_l_m,a3a_l_m,a3b_l_m,a1a_l_m,a2a_l_m,a1b_l_m,a2b_l_m,a1c_l_m,a2c_l_m
! complex*16 g_G,g_H
! complex*16 g_z(2)
! g_z=c0
! g_G=c0
! g_H=c0
! !     write(*,*)"Edmonds",lap,lbp,la,lb

! If(la.ne.lap)Then
! goto 21
! End If

! Do map=max(-lap,Mp-lbp),Min(lap,Mp+lbp)
!   mbp=Mp-map
!   call clebsch(lap*1.d0,lbp*1.d0,Lp*1.d0,map*1.d0,mbp*1.d0,Mp*1.d0,CG1)

!   Do ma=max(-la,M-lb),Min(la,M+lb)
!     mb=M-ma
!     call clebsch(la*1.d0,lb*1.d0,L*1.d0,ma*1.d0,mb*1.d0,M*1.d0,CG2)

!     If(ma.ne.map)then
!     goto 19
!     End IF
    
!     g_G=c0
! !seno*partial( Ylm  )
!     Suma1=1.d0*((2*lb+1)*(2*lb+3))
!     Suma1=1.d0*(((lb-mb+1)*(lb+mb+1))/Suma1)
!     Suma1=lb*(Suma1**0.5d0)
!     If((lbp.eq.(lb+1)).and.(mbp.eq.mb))Then
!       g_G=g_G+c1*Suma1
!     End IF

!     Suma1=1.d0*((2*lb+1)*(2*lb-1))
!     Suma1=1.d0*(((lb-mb)*(lb+mb))/Suma1)
!     Suma1=(lb+1)*(Suma1**0.5d0)
!     If((lbp.eq.(lb-1)).and.(mbp.eq.mb))Then
!       g_G=g_G-c1*Suma1
!     End IF
    
!     g_z(1)=g_z(1)+(g_G*CG1*CG2)

!     g_H=c0
! !coseno* Ylm
!     Suma1=1.d0*((2*lb+1)*(2*lb+3))
!     Suma1=1.d0*(((lb-mb+1)*(lb+mb+1))/Suma1)
!     Suma1=(Suma1**0.5d0)
!     If((lbp.eq.(lb+1)).and.(mbp.eq.mb))Then
!       g_H=g_H+c1*Suma1
!     End IF
    
!     Suma1=1.d0*((2*lb+1)*(2*lb-1))
!     Suma1=1.d0*(((lb-mb)*(lb+mb))/Suma1)
!     Suma1=(Suma1**0.5d0)
!     If((lbp.eq.(lb-1)).and.(mbp.eq.mb))Then
!       g_H=g_H+c1*Suma1
!     End IF

!     g_z(2)=g_z(2)+(g_H*CG1*CG2)
    
!     19  continue
!   End Do
! End Do

! 21 continue
! end subroutine GradienteEd_z2



! subroutine GradienteEd_1D(lap,map,la,ma,g_z)
! ! empleando las fórmulas de la página 85 del libro de Edmonds Corregido por que la deficinión de los esféricos es diferente
! use pi_module
! use complex_module
! implicit none
! integer lap,la
! integer ma,map
! real*8 Pno,Suma1,Suma2
! complex*16 g_G,g_H
! complex*16 g_z(2)
! g_z=c0
! !     write(*,*)"Edmonds",lap,lbp,la,lb


!     g_G=c0
! !seno*partial( Ylm  )
! Suma1=1.d0*((2*la+1)*(2*la+3))
! Suma1=1.d0*(((la-ma+1)*(la+ma+1))/Suma1)
! Suma1=la*(Suma1**0.5d0)
! If((lap.eq.(la+1)).and.(map.eq.ma))Then
!   g_G=g_G+c1*Suma1
! End IF

! Suma1=1.d0*((2*la+1)*(2*la-1))
! Suma1=1.d0*(((la-ma)*(la+ma))/Suma1)
! Suma1=(la+1)*(Suma1**0.5d0)
! If((lap.eq.(la-1)).and.(map.eq.ma))Then
!   g_G=g_G-c1*Suma1
! End IF

! g_z(1)=g_z(1)+(g_G)

    
! !coseno* Ylm
! g_H=c0

! Suma1=1.d0*((2*la+1)*(2*la+3))
! Suma1=1.d0*(((la-ma+1)*(la+ma+1))/Suma1)
! Suma1=(Suma1**0.5d0)
! If((lap.eq.(la+1)).and.(map.eq.ma))Then
!   g_H=g_H+c1*Suma1
! End IF

! Suma1=1.d0*((2*la+1)*(2*la-1))
! Suma1=1.d0*(((la-ma)*(la+ma))/Suma1)
! Suma1=(Suma1**0.5d0)
! If((lap.eq.(la-1)).and.(map.eq.ma))Then
!   g_H=g_H+c1*Suma1
! End IF

! g_z(2)=g_z(2)+(g_H)
    


! 21 continue
! end subroutine GradienteEd_1D



! !coeficientes de expansión!coeficientes de expansión!coeficientes de expansión!
! !coeficientes de expansión!coeficientes de expansión!coeficientes de expansión!
! !coeficientes de expansión!coeficientes de expansión!coeficientes de expansión!
! !coeficientes de expansión!coeficientes de expansión!coeficientes de expansión!
! complex*16 Function a1_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! a1_l_m=c0
!   If(l.eq.1)Then
!     If(m.eq.1)Then
!     a1_l_m=-c1*sqrt((2.d0*pi)/3.d0)
!     End If
!     If(m.eq.-1)Then
!     a1_l_m=c1*sqrt((2.d0*pi)/3.d0)
!     End If
!   End If
! End
! complex*16 Function a2_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! a2_l_m=c0
!   If(l.eq.1)Then
!     If((m.eq.1).or.(m.eq.-1))Then
!     a2_l_m=ci*sqrt((2.d0*pi)/3.d0)
!     End If
!   End If
! End
! complex*16 Function a3_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! a3_l_m=c0
!   If(l.eq.1)Then
!     If(m.eq.0)Then
!     a3_l_m=c1*2.d0*sqrt(pi/3.d0)
!     End If
!   End If
! End
! complex*16 Function a3a_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! complex*16, external :: a3_l_m
! a3a_l_m=a3_l_m(l,m)
! End
! complex*16 Function a3b_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! a3b_l_m=c0
!   If(l.eq.1)Then
!     If(m.eq.-1)Then
!     a3b_l_m=c1*2.d0*sqrt((2.d0*pi)/3.d0)
!     End If
!   End If
! End
! complex*16 Function a1a_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! a1a_l_m=c0
!   If(l.eq.1)Then
!     If(m.eq.1)Then
!     a1a_l_m=-c1*sqrt(pi/6.d0)
!     End If
!     If(m.eq.-1)Then
!     a1a_l_m=c1*sqrt(pi/6.d0)
!     End If
!   Else
!     If(m.eq.1)Then
!     a1a_l_m=-c1*sqrt((pi*(2*l+1))/(l*(l+1)))
!     End If
!     If(m.eq.-1)Then
!     a1a_l_m=c1*sqrt((pi*(2*l+1))/(l*(l+1)))
!     End If
!   End If
! a1a_l_m=a1a_l_m*(1.d0-((-1.d0)**l))*0.5d0
! End
! complex*16 Function a2a_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! a2a_l_m=c0
!   If(l.eq.1)Then
!     If((m.eq.1).or.(m.eq.-1))Then
!     a2a_l_m=ci*sqrt(pi/6.d0)
!     End If
!   Else
!     If((m.eq.1).or.(m.eq.-1))Then
!     a2a_l_m=ci*sqrt((pi*(2*l+1))/(l*(l+1)))
!     End If
!   End If
! a2a_l_m=a2a_l_m*(1.d0-((-1.d0)**l))*0.5d0
! End
! complex*16 Function a1b_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! Complex*16, external :: cgamma
! a1b_l_m=c0
!   If(l.eq.1)Then
!     If(m.eq.0)Then
!     a1b_l_m=c1*sqrt(pi/3.d0)
!     End If
!   Else
!     If(m.eq.-2)Then
!     a1b_l_m=c1*sqrt(4.d0*pi*(2*l+1))*((cgamma(c1*(l-1))/cgamma(c1*(l+3)))**0.5d0)
!     End If
!   End If
! a1b_l_m=a1b_l_m*(1.d0-((-1.d0)**l))*0.5d0
! End

! complex*16 Function a2b_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! Complex*16, external :: cgamma
! a2b_l_m=c0
!   If(l.eq.1)Then
!     If(m.eq.0)Then
!     a2b_l_m=-ci*sqrt(pi/3.d0)
!     End If
!   Else
!     If(m.eq.-2)Then
!     a2b_l_m=ci*sqrt(4.d0*pi*(2*l+1))*((cgamma(c1*(l-1))/cgamma(c1*(l+3)))**0.5d0)
!     End If
!   End If
! a2b_l_m=a2b_l_m*(1.d0-((-1.d0)**l))*0.5d0
! End

! complex*16 Function a1c_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! a1c_l_m=c0
! !   If(l.eq.1)Then
!     If((m.eq.1).or.(m.eq.-1))Then
!     a1c_l_m=ci*(1.d0-((-1.d0)**l))*0.5d0*sqrt((pi*(2*l+1))/(1.d0*l*(l+1)))
!     End If
! !   End If
! End

! complex*16 Function a2c_l_m(l,m)
! use pi_module
! use complex_module
! implicit none
! integer l,m
! a2c_l_m=c0
! !   If(l.eq.1)Then
!     If(m.eq.1)Then
!     a2c_l_m=-c1*(1.d0-((-1.d0)**l))*0.5d0*sqrt((pi*(2*l+1))/(1.d0*l*(l+1)))
!     ElseIf(m.eq.-1)Then
!     a2c_l_m=c1*(1.d0-((-1.d0)**l))*0.5d0*sqrt((pi*(2*l+1))/(1.d0*l*(l+1)))
!     End If
! !   End If
! End

! !coeficientes de expansión!coeficientes de expansión!coeficientes de expansión!
! !coeficientes de expansión!coeficientes de expansión!coeficientes de expansión!
! !coeficientes de expansión!coeficientes de expansión!coeficientes de expansión!
! !coeficientes de expansión!coeficientes de expansión!coeficientes de expansión!


! ! da la parte angular del elemento de matriz de ri vector.
! subroutine vector_r1(Lp,Mp,lap,lbp,L,M,la,lb,v_)
! use pi_module
! use complex_module
! implicit none
! integer Lp,Mp,lap,lbp,L,M,la,lb,le
! integer ma,mb,map,mbp,me
! real*8 CG1,CG2,Pno
! real*8, external :: I3
! complex*16, external :: a1_l_m,a2_l_m,a3_l_m,a3a_l_m,a3b_l_m,a1a_l_m,a2a_l_m,a1b_l_m,a2b_l_m,a1c_l_m,a2c_l_m
! complex*16 v_x,v_y,v_z
! complex*16 v_(3),Suma1,Suma2
! ! write(*,*) "vector_r1"
! v_=c0
! v_x=c0
! v_y=c0
! v_z=c0

! If(lb.ne.lbp)Then
! ! write(*,*) "lb.ne.lbp"
! goto 21
! End If
! ! write(*,*) max(-lap,Mp-lbp),Min(lap,Mp+lbp)
! Do map=max(-lap,Mp-lbp),Min(lap,Mp+lbp)
!   mbp=Mp-map
!   call clebsch(lap*1.d0,lbp*1.d0,Lp*1.d0,map*1.d0,mbp*1.d0,Mp*1.d0,CG1)
! ! write(*,*) CG1,max(-la,M-lb),Min(la,M+lb),M-Mp+map
!   Do ma=max(-la,M-lb),Min(la,M+lb)
! !     If(ma.ne.(M-Mp+map))Then
! !       goto 19
! !     End If

!     mb=M-ma
    
    
! !     If(mb.ne.mbp)Then
! !       goto 19
! !     End If
    
    
    
    
!     call clebsch(la*1.d0,lb*1.d0,L*1.d0,ma*1.d0,mb*1.d0,M*1.d0,CG2)
    
! !     write(*,*) map,ma,mbp,mb,CG1*CG2
    
    
! ! write(*,*) CG2,int(abs(la-lap)*1.001d0),(la+lap)
!     Do le=int(abs(la-lap)*1.001d0),(la+lap)
!       Do me=-le,le         !map-ma-1,map-ma			!


! ! 	If(CG1*CG2.eq.0.d0)Then
! ! 	  goto	20
! ! 	End If

! 	If(mbp.eq.mb)Then
! 	  Suma1=CG1*CG2*I3(lap,map,la,ma,le,me)

! 	  v_x=v_x+c1*a1_l_m(le,me)*Suma1
! 	  v_y=v_y+c1*a2_l_m(le,me)*Suma1
! 	  v_z=v_z+c1*a3_l_m(le,me)*Suma1
! ! write(*,*)Suma1,v_z
! 	End If

! 	20 continue


!       End Do
!     End Do

!     19 continue

!   End Do
! End Do

! v_(1)=v_x
! v_(2)=v_y
! v_(3)=v_z

! 21 continue
! end subroutine vector_r1

! subroutine vector_r2(Lp,Mp,lap,lbp,L,M,la,lb,v_)
! use pi_module
! use complex_module
! implicit none
! integer Lp,Mp,lap,lbp,L,M,la,lb,le
! integer ma,mb,map,mbp,me
! real*8 CG1,CG2,Pno
! real*8, external :: I3
! complex*16, external :: a1_l_m,a2_l_m,a3_l_m,a3a_l_m,a3b_l_m,a1a_l_m,a2a_l_m,a1b_l_m,a2b_l_m,a1c_l_m,a2c_l_m
! complex*16 v_x,v_y,v_z
! complex*16 v_(3),Suma1,Suma2
! ! write(*,*) "vector_r2"
! v_=c0
! v_x=c0
! v_y=c0
! v_z=c0

! If(la.ne.lap)Then
! ! write(*,*) "la.ne.lap"

! goto 21
! End If
! ! write(*,*) max(-lap,Mp-lbp),Min(lap,Mp+lbp)
! Do map=max(-lap,Mp-lbp),Min(lap,Mp+lbp)
!   mbp=Mp-map
!   call clebsch(lap*1.d0,lbp*1.d0,Lp*1.d0,map*1.d0,mbp*1.d0,Mp*1.d0,CG1)
! ! write(*,*)CG1,max(-la,M-lb),Min(la,M+lb)
!   Do ma=max(-la,M-lb),Min(la,M+lb)
! !     If(ma.ne.(M-Mp+map))Then
! !       goto 19
! !     End If

!     mb=M-ma
    
! !     If(ma.ne.map)Then
! !       goto 19
! !     End If
    
!     call clebsch(la*1.d0,lb*1.d0,L*1.d0,ma*1.d0,mb*1.d0,M*1.d0,CG2)
! ! write(*,*) CG2,int(abs(lb-lbp)*1.001d0),(lb+lbp)
!     Do le=int(abs(lb-lbp)*1.001d0),(lb+lbp)
!       Do me=-le,le         !map-ma-1,map-ma			!


! ! 	If(CG1*CG2.eq.0.d0)Then
! ! 	  goto	20
! ! 	End If

! 	If(map.eq.ma)Then
! 	  Suma1=CG1*CG2*c1*I3(lbp,mbp,lb,mb,le,me)

! 	  v_x=v_x+c1*a1_l_m(le,me)*Suma1
! 	  v_y=v_y+c1*a2_l_m(le,me)*Suma1
! 	  v_z=v_z+c1*a3_l_m(le,me)*Suma1
! ! write(*,*)Suma1,v_z
! 	End If

! 	20 continue


!       End Do
!     End Do

!     19 continue

!   End Do
! End Do

! v_(1)=v_x
! v_(2)=v_y
! v_(3)=v_z

! 21 continue
! end subroutine vector_r2

!     real*8 function A_exp(l)
!     implicit none
!     integer*4 l
! 	real*8 a(0:50)
! 	data a /0.d0,3.4098905781729165d0,0.d0,0.7974151418694725d0,0.d0,&
!        0.39513164604699597d0,0.d0,0.24625469619275697d0,0.d0,0.1721626796366226d0,0.d0,&
!        0.12903703307480943d0,0.d0,0.10134687461820338d0,0.d0,0.08232391039179625d0,0.d0,&
!        0.06859173160629006d0,0.d0,0.05829603062216876d0,0.d0,0.05034217078448181d0,0.d0,&
!        0.04404648909982984d0,0.d0,0.03896220899510459d0,0.d0,0.0347862811248568d0,0.d0,&
!        0.03130663651362859d0,0.d0,0.028370887962959113d0,0.d0,0.02586700823335441d0,0.d0,&
!        0.02371099050668935d0,0.d0,0.02183873292840546d0,0.d0,0.020200562156658462d0,0.d0,&
!        0.018757452859631544d0,0.d0,0.017478364689172736d0,0.d0,0.016338332160540155d0,0.d0,&
!        0.015317072051757774d0,0.d0,0.014397953004342498d0,0.d0/
! 	if(l.gt.50) then
! 	    write(*,*) 'too large l in A_exp'
! 	    return
! 	end if
! 	A_exp = a(l)
! 	return
! 	end

!     real*8 function A_sinm1(l)
!     implicit none
!     integer*4 l
! 	real*8 a(0:50)
! 	data a /5.568327996831707d0,0.d0,3.1127899804827326d0,0.d0,&
!         2.349138373663376d0,0.d0,1.9606339952322658d0,0.d0,1.7165854186432425d0,0.d0,&
!         1.5453812002343439d0,0.d0,1.4168336402842365d0,0.d0,1.3157654786435473d0,0.d0,&
!         1.2336132106222124d0,0.d0,1.1651336257776372d0,0.d0,1.106914338348481d0,0.d0,&
!         1.0566266725643294d0,0.d0,1.0126200942089252d0,0.d0,0.9736878647948234d0,0.d0,&
!         0.9389245961788506d0,0.d0,0.9076359558645707d0,0.d0,0.879279370448314d0,0.d0,&
!         0.8534238913292798d0,0.d0,0.8297223113464376d0,0.d0,0.8078913457175132d0,0.d0,&
!         0.7876972588784051d0,0.d0,0.7689452530268338d0,0.d0,0.751471507813666d0,0.d0,&
!         0.7351371224677103d0,0.d0,0.7198234454392937d0,0.d0,0.7054284310176197d0/
! 	if(l.gt.50) then
! 	    write(*,*) 'too large l in A_sinm1'
! 	    return
! 	end if
! 	A_sinm1 = a(l)
! 	return
! 	end



!     real*8 function A_ctag(l)
!     implicit none
!     integer*4 l
! 	real*8 a(0:50)
! 	data a /0.d0,4.3416075273496055d0,0.d0,2.7074679791968332d0,0.d0,2.1465482117262193d0,&
!        0.d0,1.8346640370504368d0,0.d0,1.6287726861141123d0,0.d0,1.4797279749246832d0,0.d0,&
!        1.3653724013755d0,0.d0,1.274031476987412d0,0.d0,1.1988869722499382d0,0.d0,&
!        1.1356524370552348d0,0.d0,1.0814792022839936d0,0.d0,1.0343900209075132d0,0.d0,&
!        0.9929636303916116d0,0.d0,0.956148568151405d0,0.d0,0.9231479540743802d0,0.d0,&
!        0.8933453283642951d0,0.d0,0.8662553002680455d0,0.d0,0.8414897568058258d0,0.d0,&
!        0.8187341459605134d0,0.d0,0.7977304660700844d0,0.d0,0.7782648293444512d0,0.d0,&
!        0.7601582132836867d0,0.d0,0.7432594770310473d0,0.d0,0.7274400149419237d0,0.d0,&
!        0.7125896122154537d0,0.d0/
! 	if(l.gt.50) then
! 	    write(*,*) 'too large l in A_ctag'
! 	    return
! 	end if
! 	A_ctag = a(l)
! 	return
! 	end
! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''





!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!
!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!
!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!
!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!ENERGIA NO ORTOGONAL!


!!!VECTOR C2!!!VECTOR C2!!!VECTOR C2!!!VECTOR C2!!!VECTOR C2!!!VECTOR C2!!!VECTOR C2



!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS



!	ARTURO QUIRANTES SIERRA
!	Department of Applied Physics, Faculty of Sciences
!	University of Granada, 18071 Granada (SPAIN)
!	http://www.ugr.es/local/aquiran/codigos.htm
!	aquiran@ugr.es
!	Last update: 20 May 2.003
!	Subroutine NED
!	to calculate Clebsch-Gordan coefficients
!	You need to add a "NED(AJ,BJ,CJ,AM,BM,CM,CG)" in your main routine
!	Input:
!		AJ,BJ,CJ,AM,BM,CM (the usual Clebsch-Gordan indices)
!	Output:
!		CG=C-G(AJ,BJ,CJ,AM,BM,CM)
                        !  (J,J1,J2,M,M1,M2,CG)

	SUBROUTINE clebsch(AJ,BJ,CJ,AM,BM,CM,CG)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION Q(500,500)
	DOUBLE PRECISION CG
	DOUBLE PRECISION AJ,BJ,CJ,AM,BM,CM
	INTEGER ZZ
	ZZ=MAX(2*AJ+1,2*BJ+1,2*CJ+1,AJ+BJ+CJ,AJ+AM,BJ+BM,CJ+CM)+2
	DO 2 I=1,ZZ
		Q(I,1)=1.D0
		Q(I,I)=1.D0
2	CONTINUE
	DO 3 I=2,ZZ-1
	DO 3 K=2,I
		Q(I+1,K)=Q(I,K-1)+Q(I,K)
3	CONTINUE
	CG=0.D0
	JA=AJ+AM+1.01D0
	MA=AJ-AM+1.01D0
	JB=BJ+BM+1.01D0
	MB=BJ-BM+1.01D0
	JC=CJ+CM+1.01D0
	MC=CJ-CM+1.01D0
	LA=BJ+CJ-AJ+1.01D0
	LB=CJ+AJ-BJ+1.01D0
	LC=AJ+BJ-CJ+1.01D0
	LT=AJ+BJ+CJ+1.01D0
	D=DABS(AM+BM-CM)-0.01D0
	IF (D) 10,10,20
10	LD=MIN0(JA,JB,JC,MA,MB,MC,LA,LB,LC)
	IF (LD) 20,20,30
30	JA2=AJ+AJ+AM+AM
	JB2=BJ+BJ+BM+BM
	JC2=CJ+CJ-CM-CM
	I2=JA2+JB2+JC2-JA2/2*2-JB2/2*2-JC2/2*2
	IF (I2) 20,40,20
40	FN=Q(JA+MA-1,LC)/Q(LT,JC+MC-1)
	FN=FN*Q(JB+MB-1,LC)/Q(LT+1,2)
	FN=FN/Q(JA+MA-1,JA)
	FN=FN/Q(JB+MB-1,JB)
	FN=FN/Q(JC+MC-1,JC)
	K0=MAX(0,LC-JA,LC-MB)+1
	K1=MIN(LC,MA,JB)
	X=0.D0
	DO 50 K=K0,K1
		X=-X-Q(LC,K)*Q(LB,MA-K+1)*Q(LA,JB-K+1)
50	CONTINUE
	IP=K1+LB+JC
	P=1-2*(IP-IP/2*2)
	CG=P*X*DSQRT(FN)
!	What weŽve calculated is a Wigner 3-j coefficient
!	Next, weŽll turn it into a Clebsch-Gordan coefficient
	CG=CG*DSQRT(2*CJ+1)*(-1)**IDNINT(AJ-BJ-CM)
20	CONTINUE
	RETURN
	END


	SUBROUTINE APU(IAJ,IBJ,ICJ,IAM1,IAM2,IBM,CGN)
	double precision cc1,cc2,dd1,dd2,cgn(500),CX(500),C1,SC,CMAX,CMAX2
	INTEGER IAJ,IBJ,ICJ,IBM,IAM1,IAM2,IAM,IAMEDIO
	INTEGER DN1,KP,KP2,IAMA,IAMB,DNTOPE,IAMIN,IAMAX
!	First, let's check that icj and ibm falls within limits
	if(icj.lt.abs(iaj-ibj).or.icj.gt.(iaj+ibj)) go to 100
	if(abs(ibm).gt.ibj) go to 200
!	IAMIN, IAMAX are the lowest and highest possible values for am
	IAMIN=-MIN(IAJ,ICJ+IBM)
	IAMAX=MIN(IAJ,ICJ-IBM)
!	IAMIN, IAMAX might be such that CG=0
!	If so happens, let's jump up till CG<>0
	KP=MAX(0,IAMIN-IAM1)
	KP2=MAX(0,IAM2-IAMAX)
	DNTOPE=10
!	Calculation of every possible CG value
	IAMA=IAMIN
	IAMB=IAMAX
	DN1=IAMAX-IAMIN
	CX(1)=1.0D0
	CMAX=CX(1)
	SC=CX(1)*CX(1)
	if(dn1.eq.0) go	to 100
	dd2=dfloat(iaj*(iaj+1)-ibj*(ibj+1)+icj*(icj+1))
	cc2=dfloat((iaj-iama+1)*(iaj+iama))
	cc2=cc2*dfloat((icj-iama-ibm+1)*(icj+iama+ibm))
	cc2=-sqrt(cc2)
	IF(DN1.GE.DNTOPE) IAMB=IAMIN+INT(DN1/2)-1
!	If there are not many CG to calculate, upward recurrence alone will do
	iam=iama
5	iam=iam+1
		icm=iam+ibm
		if(icm.gt.icj) go to 8
		cc1=cc2
		cc2=dfloat((iaj-iam+1)*(iaj+iam))
		cc2=cc2*dfloat((icj-icm+1)*(icj+icm))
		cc2=-sqrt(cc2)
		if(-icm.gt.icj) go to 8
		dd1=dd2-2.0d0*dfloat((iam-1)*(icm-1))
		CX(IAM-IAMIN+1)=-CX(IAM-IAMIN)*DD1/CC2
		if(dn1.eq.1) then
			SC=SC+CX(IAM-IAMIN+1)*CX(IAM-IAMIN+1)
			go to 100
		end if
		if(iam.eq.(iama+1)) go to 7
		CX(IAM-IAMIN+1)=CX(IAM-IAMIN+1)-CX(IAM-IAMIN-1)*CC1/CC2
7		SC=SC+CX(IAM-IAMIN+1)*CX(IAM-IAMIN+1)
		CMAX2=ABS(CX(IAM-IAMIN+1))
		IF(CMAX2.GT.CMAX) CMAX=CMAX2
!	Now let's see if we need go on with recurrence
		if(iam.lt.iamb) go to 5
!	CX(middle) is used for re-scaling.  If equals zero, then disaster!
!	So let us just check and, if zero, go for the next one.
!	Due to round-off errors, we have to be a bit more tolerant: X<10^-14
!	(X being the smallest/largest C-G ratio in our series)
!	The problem is, for high aj, bj... values, CG can have very small values
!	But then you would need extended precision anyway, so you might change
!	the tolerance criterion to, say, (c2max/cmax)<10¯30
	if(dn1.ge.dntope.and.iam.ge.iamb.and.(cmax2/cmax).le.(1.0D-14)) go to 5
8	continue
	if(dn1.lt.dntope) go to 100
	iamedio=iam
	C1=CX(IAMEDIO-IAMIN+1)
	SC=SC-C1*C1
!	Downward recurrence
	IAMA=IAMEDIO
	iamb=iamax
	CX(IAMB-IAMIN+1)=1.0D0
	cc1=dfloat((iaj-iamb)*(iaj+iamb+1))
	cc1=cc1*dfloat((icj-iamb-ibm)*(icj+iamb+ibm+1))
	cc1=-sqrt(cc1)
	do iam=iamb-1,iama,-1
		icm=iam+ibm
		if(icm.gt.icj) go to 100
		cc2=cc1
		cc1=dfloat((iaj-iam)*(iaj+iam+1))
		cc1=cc1*dfloat((icj-icm)*(icj+icm+1))
		cc1=-dsqrt(cc1)		
		if(-icm.gt.icj) go to 18
		dd1=dd2-2.0d0*dfloat((iam+1)*(icm+1))
		CX(IAM-IAMIN+1)=-CX(IAM-IAMIN+2)*DD1/CC1
		if(iam.eq.(iamb-1)) go to 18
		CX(IAM-IAMIN+1)=CX(IAM-IAMIN+1)-CX(IAM-IAMIN+3)*CC2/CC1
18	end do
	C1=C1/CX(IAMEDIO-IAMIN+1)
!	Now, rescaling
	do iam=IAMEDIO,IAMAX
		cx(iam-iamin+1)=cx(iam-iamin+1)*c1
		sc=sc+cx(iam-iamin+1)*cx(iam-iamin+1)
	end do
!	Here's where the combined (upward+downward) recurrence ends
100	CONTINUE
	SC=(dfloat(2*ICJ+1)/dfloat(2*IBJ+1))/SC
        SC=SQRT(SC)
!	Now, let's set the sign for SC
	if(cx(iamax-iamin+1).ne.abs(cx(iamax-iamin+1))) SC=-SC
	if(mod((iaj+iamax),2).ne.0) SC=-SC
	if(kp.ne.0) then
		do iam=0,kp-1
			cgn(iam+1)=0.0d0				
		end do
	end if
	if(kp2.ne.0) then
		do iam=iam2-kp2+1,iam2
			cgn(iam-iam1+1)=0.0d0
		end do
	end if
	do iam=iam1+kp,iam2-kp2
		cgn(iam-iam1+1)=cx(iam-iamin+1)*SC
	end do
200	CONTINUE
	RETURN
	END

!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS!FUNCIONES Y ACCESORIOS
