
!real*8 Function HarmonicOscillator(n, ome, x)
real*8 Function HO(n, ome, x)
implicit none
integer, intent(in) :: n
real*8 , intent(in) :: ome, x
REAL*8, external :: FACT
real*8 p(1,0:n)
real*8 x0,xm(1)

x0 = x*sqrt(ome)
xm = x0
call h_polynomial_value ( 1, n, xm, p )

HO = p(1,n)*exp(-x0*x0*0.5d0)*((ome/3.14159265d0)**0.25d0)
HO = HO/sqrt((2**n)*fact(n))

End Function HO



! real*8 Function FreeCoulomb(n, L, r)
! use ModuleParameters
! implicit none
! real*8  :: r
! integer :: L, n
! real*8  :: f, fp, g, gp, En
! En = Energias(n,L)
! call scoul(Ztilde*sqrt(me),En,0,(Rmax-r)*sqrt(me),f,fp,g,gp,err)
! FreeCoulomb=f*NormCoul(n,L)
! End Function FreeCoulomb





! program example
!     use root_finder
!     implicit none
!     real(dp) :: root
!
!     ! Encuentra el primer cero de J0(x) en [2.0, 5.0]
!     root = find_root_brent(bessel_j0, 2.0_dp, 5.0_dp)
!     print *, "Primer cero de J0(x) ≈ ", root
!
! contains
!
!     function bessel_j0(x) result(j0)
!         real(dp), intent(in) :: x
!         real(dp) :: j0
!         ! Usa una implementación de J0(x) (puede ser intrínseca, GSL, o SLATEC)
!         j0 = BESSEL_J0(x)  ! Si usas Intel/gfortran moderno
!         ! O alternativamente: call dbesj(x, 0.0_dp, 1, j0)  ! SLATEC
!     end function bessel_j0
!
! end program example



!!!indice a partir de pasarle una tupla (n1 n2 n3 n4...)
    integer function Index_photonic(tupla, nmax) result(indice)
        implicit none
        integer, intent(in) :: tupla(:)
        integer, intent(in) :: nmax
        integer :: m, i, k
        indice = 0
        m = size(tupla)
        do k = 1, m
            indice = indice + tupla(k) * (nmax+1)**(k-1)
        end do
    end function Index_photonic

    function fock_photonic(indice, nmax, vector_size) result(tupla)  !!hay que pasarle el indice, el nro max de fotones y la cantidad de modos
        implicit none
        integer, intent(in) :: indice, nmax, vector_size
        integer :: tupla(vector_size)
        integer :: b, N, r, j, digit_count
        
        ! Inicializar
        tupla = 0
        N = indice
        b = nmax + 1
        digit_count = vector_size
        
        ! Conversión CORRECTA a base b
        do j = vector_size, 1, -1
            r = mod(N, b)
            tupla(j) = r
            N = N / b
            if (N == 0) exit
        end do
        
    end function fock_photonic

!!funciona, da bien ambos numeros





   








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! real*8 Function SphericalBesselJ(L, r)
! use ModuleConstants
! use ModuleParameters
! implicit none
! real*8        :: r
! integer       :: L
! real*8        :: f
! integer, save :: jo
! ! write(*,*) "aca",rsbound,Nrbound,
! Call hunt(rsbound,Nrbound,r,jo)
! ! write(*,*) "aca 2"
! !call locate(rad,Nrbound,r,jo) para evaluación random
! f=jsa(L,jo)+r*(jsb(L,jo)+r*(jsc(L,jo)+r*jsd(L,jo)))
! ! write(*,*) "aca 3"
! SphericalBesselJ=f
! End Function SphericalBesselJ
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine WriteBasis
! use ModuleConstants
! use ModuleParameters
! implicit none
! integer           :: L, n
! Real*8            :: r, dr, func, func2
! character(len=8)  :: fmt1, fmt2 ! format descriptor
! character(5)      :: filenumber1 , filenumber2
! real*8,  external :: BoundCoulomb, FreeCoulomb
! logical           :: dir_e
! inquire(file='./basis/test', exist=dir_e)
! if (dir_e) then
! !   write(*,*) "dir exists!",dir_e
! else
!   ! workaround: it calls an extern program...
! !   write(*,*) "dir does not exists!",dir_e
!   call system('mkdir -p ./basis')
!   Open(unit=10,file="basis/test")
!   close(10)
! end if
! fmt1 = '(I3.3)'
! fmt2 = '(I2.2)'
! dr = Rmax/10000.d0
! Do L = 0, Lmax
!    write (filenumber2,fmt2) L
!    Do n = 1, Nmax
!       write (filenumber1,fmt1) n
!       func2 = 0.d0
!       open(2, file = 'basis/CoulombBox_L'//trim(filenumber2)//'_n'//trim(filenumber1)//'.dat')
!       do r=0.d0,Rmax,dr
!          func = FreeCoulomb(n, L, r)
!          func2 = func2+func*func*dr
!          write(2,*) r,func
!       enddo
!       write(*,*) "check norma free", L,n,func2
!       close(2)
!    enddo
! enddo
! Do L = 0, Lmax
!     write (filenumber2,fmt2) L
!    Do n = L+1,Nbounds
!       func2 = 0.d0
!       write (filenumber1,fmt1) n
!       open(2, file = 'basis/BoundBox_L'//trim(filenumber2)//'_n'//trim(filenumber1)//'.dat')
!       do r=0.d0,Rmax,dr
!          func = BoundCoulomb(n, L, r)
!          func2 = func2+func*func*dr
!          write(2,*) r,func
!       enddo
!       write(*,*) "check norma bound", L,n,func2
!       close(2)
!    enddo
! enddo
! endsubroutine
!
! subroutine WriteBessels
! use ModuleConstants
! use ModuleParameters
! implicit none
! integer           :: L
! Real*8            :: r, dr, func
! character(len=8)  :: fmt2 ! format descriptor
! character(5)      :: filenumber2
! real*8,  external :: SphericalBesselJ
! logical           :: dir_e
! inquire(file='./bessels/test', exist=dir_e)
! if (dir_e) then
! !   write(*,*) "dir exists!",dir_e
! else
!   ! workaround: it calls an extern program...
! !   write(*,*) "dir does not exists!",dir_e
!   call system('mkdir -p ./bessels')
!   Open(unit=10,file="bessels/test")
!   close(10)
! end if
! fmt2 = '(I2.2)'
! dr = Rmax/10000.d0
! Do L = 0, LmaxJ
!    write (filenumber2,fmt2) L
!    open(2, file = 'bessels/SphBessel_L'//trim(filenumber2)//'.dat')
!    do r=0.d0,Rmax,dr
!       func = SphericalBesselJ(L, r)
!       write(2,*) r,func
!    enddo
!    close(2)
! enddo
!
! endsubroutine
!
!
!
!
REAL*8 FUNCTION LegendreP(l,m,x)	!plgndr(l,m,x)
implicit none
INTEGER l,m
double precision x
!Computes the associated Legendre polynomial Plm (x). Here m and l are integers satisfying
!0 ≤ m ≤ l, while x lies in the range −1 ≤ x ≤ 1.
INTEGER i,ll
REAL*8 fact,pll,pmm,pmmp1,somx2
!if (((m .lt. 0) .or. (m .gt. l)) .or. (abs(x).gt.1.d0)) write(*,*) "’bad arguments in plgndr: Ver Polinomios de Legendre’"
! write(*,*) "Legendre l,m,x",l,m,x
! if (((m .lt. 0) .or. (m .gt. l)) .or. (abs(x).gt.1.d0))Then
! write(*,*) "Argumento no válido en los polinomios de Legendre",l,m,abs(x)
! stop
! end if
pmm=1.d0			!Compute P_{m}^{m} .
if(m.gt.0) then
    somx2=sqrt((1.d0-x)*(1.d0+x))
    fact=1.d0
    do i=1,m
         pmm=-pmm*fact*somx2
         fact=fact+2.d0
    enddo
endif
if(l .eq. m) then
    LegendreP=pmm
else
pmmp1=x*(2.d0*m+1.d0)*pmm		!Compute P_{m+1}^{m}
    if(l .eq. (m+1)) then
         LegendreP=pmmp1
    else		!Compute P_{l}^{m}, l > m + 1.
         do ll=m+2,l
            pll=(x*(2.d0*ll-1.d0)*pmmp1-(ll+m-1.d0)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
        enddo
        LegendreP=pll
    endif
endif
return
END FUNCTION LegendreP
!
!
! complex*16 FUNCTION factGammaCoc(l,m)
! use ModuleConstants
! implicit none
! integer l,m
! complex*16, external :: cgamma,cgammar
! !c1*(sqrt(factGamma(l-m))/sqrt(factGamma(l+m)))
! factGammaCoc=c1*sqrt(real(cgamma(c1*(l-m+1))*cgammar(c1*(l+m+1)),kind=8))
! End Function
!
complex*16 FUNCTION SphericalH(l,m,theta,varphi)
use complex_module
use pi_module
implicit none
INTEGER l,m,k,kk
double precision theta,varphi
complex*16, external :: cgamma,cgammar
real*8, external :: LegendreP,fact
complex*16, external :: factGammaCoc


SphericalH=c1*(sqrt(2.d0*l+1.d0)/sqrt(4.d0*pi))*((cgamma((l-m+1)*c1)/cgamma((l+m+1)*c1))**0.5d0)*exp(ci*m*varphi)
If (m .ge. 0) Then


    SphericalH=SphericalH*LegendreP(l,m,1.d0*cos(theta))



Else
    SphericalH=SphericalH*((-1)**m)*LegendreP(l,-m,1.d0*cos(theta))*(cgamma((l+m+1)*c1)/cgamma((l-m+1)*c1))
End If

If(theta.lt.0.d0)Then
SphericalH=SphericalH*((-1)**m)
End If


! kk=l
! k=m
! If (m .ge. 0) Then
! ! SphericalH=c1*(sqrt((2.d0*l+1.d0)*factGamma(l-m))/sqrt(4.d0*pi*factGamma(l+m)))
!     SphericalH=c1*(sqrt(2.d0*l+1.d0)/sqrt(4.d0*pi))*factGammaCoc(l,m)
!     SphericalH=SphericalH*LegendreP(kk,k,1.d0*cos(theta))*exp(ci*m*varphi)
! !     write(*,*) "1111l,m,theta,varphi",l,m,theta,varphi
! Else
! !     SphericalH=c1*(sqrt((2.d0*l+1.d0)*factGamma(l-(-m)))/sqrt(4.d0*pi*factGamma(l+(-m))))
!     SphericalH=c1*(sqrt(2.d0*l+1.d0)/sqrt(4.d0*pi))*factGammaCoc(l,-m)
! !     write(*,*) "SphericalH",SphericalH
!     SphericalH=SphericalH*LegendreP(kk,-k,1.d0*cos(theta))*exp(ci*(-m)*varphi)
!     SphericalH=((-1)**m)*(real(SphericalH)-ci*real(-ci*SphericalH))
! End If
END FUNCTION SphericalH
!
!
! real*8 Function I3(l1,m1,l2,m2,l3,m3)
! use ModuleConstants
! implicit none
! integer l1,m1,l2,m2,l3,m3
! real*8 CG1,CG2
! call clebsch(l1*1.d0,l2*1.d0,l3*1.d0,0.d0,0.d0,0.d0,CG1)
! call clebsch(l1*1.d0,l2*1.d0,l3*1.d0,-m1*1.d0,m2*1.d0,-m3*1.d0,CG2)
! I3=((-1.d0)**(m3+m1))*dsqrt(((2*l1+1)*(2*l2+1))/(4.d0*pi*(2*l3+1)))*CG1*CG2
! End
!
!
! subroutine WriteBasis2D
! use ModuleConstants
! use ModuleParameters
! implicit none
! integer             :: L, n, ic, id
! Real*8              :: rho2, z2, r
! character(len=8)    :: fmt1, fmt2 ! format descriptor
! character(5)        :: filenumber1 , filenumber2
! real*8,  external   :: BoundCoulomb, FreeCoulomb
!
! complex*16,  external :: SphericalH
! logical             :: dir_e
! complex*16, allocatable :: RePsi(:,:)
! real*8, allocatable :: rhopl(:), zpl(:)
!
! write(*,*) "WriteBasis2D"
! inquire(file='./Basis2D/test', exist=dir_e)
! if (dir_e) then
! !   write(*,*) "dir exists!",dir_e
! else
!   call system('mkdir -p ./Basis2D')
!   Open(unit=10,file="Basis2D/test")
!   close(10)
! end if
!
!
!
! Rplot  = 15.d0
! Zplot  = 30.d0
! Zshift = 10.d0
! Nplot  = 79 !El Rho lo divido en esta cantidad, tengo valores hacia arriba y abajo. +1 es la línea rho=0.d0
! NRplot = (2*Nplot+1)
! NZplot = Nplot*(Rplot/Zplot)!NRplot
!
!
!
! allocate(RePsi(NZplot,NRplot),rhopl(NRplot),zpl(NZplot))
! rhopl   = 0.d0
! zpl     = 0.d0
! Do ic = 1, NRplot
!    rhopl(ic) = -Rplot + (2.d0*Rplot*(ic-1)/(1.d0*(NRplot-1)))
! Enddo
! Do ic = 1, NZplot
!    zpl(ic)   = -Zplot + (2.d0*Zplot*(ic-1)/(1.d0*(NZplot-1))) + Zshift
! Enddo
!
! fmt1 = '(I3.3)'
! fmt2 = '(I2.2)'
!
! Do L = 0, Lmax
!    write (filenumber2,fmt2) L
!    write(*,*) "L continuos",L
!    Do n = 1, Nmax
!       write (filenumber1,fmt1) n
!       RePsi = c0
!
!          Do id = 1, NZplot
!             z2 = zpl(id)*zpl(id)
!                   Do ic = 1, NRplot
!          rho2 = rhopl(ic)*rhopl(ic)
!             r = (rho2 + z2)**0.5d0
!             RePsi(id,ic) = FreeCoulomb(n, L, r)*(r**(-1.d0))*SphericalH(L,0,ATAN2(zpl(id),rhopl(ic)),0.d0)
!          Enddo
!       Enddo
!       open(2, file = 'Basis2D/CoulombBox_L'//trim(filenumber2)//'_n'//trim(filenumber1),form='unformatted')
!       write(2) RePsi
!       close(2)
!    enddo
!    write(*,*) "L ligados",L
!    Do n = L+1,Nbounds
!       write (filenumber1,fmt1) n
!       RePsi = c0
!
!          Do id = 1, NZplot
!             z2 = zpl(id)*zpl(id)
!                   Do ic = 1, NRplot
!          rho2 = rhopl(ic)*rhopl(ic)
!             r = (rho2 + z2)**0.5d0
!             RePsi(id,ic) = BoundCoulomb(n, L, r)*(r**(-1.d0))*SphericalH(L,0,ATAN2(zpl(id),rhopl(ic)),0.d0)
!          Enddo
!       Enddo
!       open(2, file = 'Basis2D/BoundBox_L'//trim(filenumber2)//'_n'//trim(filenumber1),form='unformatted')
!       write(2) RePsi
!       close(2)
!    enddo
! enddo
!
! write(*,*) "aaa"
!
! open(2, file = 'Basis2D/rho',form='unformatted')
! write(2) rhopl
! close(2)
! open(2, file = 'Basis2D/z',form='unformatted')
! write(2) zpl
! close(2)
! write(*,*) "bbbb"
! deallocate(RePsi,rhopl,zpl)
! write(*,*) "ccc"
! open(2, file = 'Basis2D/constants',form='unformatted')
! write(2) Nplot
! write(2) NRplot
! write(2) NZplot
! write(2) Rplot
! write(2) Zplot
! write(2) Zshift
! close(2)
! write(*,*) "cccsssss"
! endsubroutine WriteBasis2D
!
!
!
!
!
!
! subroutine ComputeFunction(ti, to, dt)
! use ModuleConstants
! use ModuleParameters
! implicit none
! Real*8              :: ti, to, dt, t
! integer             :: L, n, ic, id
! Real*8              :: rho2, z2, r
! character(len=8)    :: fmt1, fmt2, fmt3 ! format descriptor
! character(5)        :: filenumber1 , filenumber2
! character(8)        :: timeCh
! real*8,  external   :: BoundCoulomb, FreeCoulomb
! logical             :: dir_e
! complex*16, allocatable :: Basis(:,:)
! real*8, allocatable :: RePsi(:,:), ImPsi(:,:)
! real*8, allocatable :: rhopl(:), zpl(:)
! inquire(file='./Solution/test', exist=dir_e)
! if (dir_e) then
! !   write(*,*) "dir exists!",dir_e
! else
!   call system('mkdir -p ./Solution')
!   Open(unit=10,file="Solution/test")
!   close(10)
! end if
!
! allocate(Cn_f(1:Nmax,0:Lmax,-Lmax:Lmax), Cn_b(1:Nbounds,0:Lmax,-Lmax:Lmax))
! open(2, file = 'coefficients/coefficientsFree',form='unformatted')
! read(2) Cn_f
! close(2)
! open(2, file = 'coefficients/coefficientsBound',form='unformatted')
! read(2) Cn_b
! close(2)
!
! fmt1 = '(I3.3)'
! fmt2 = '(I2.2)'
! fmt3 = '(F5.2)'
! write(timeCh,fmt3) t
! write(*,*) "t",t, "......",timeCh
!
! open(2, file = 'Basis2D/constants',form='unformatted')
! read(2) Nplot
! read(2) NRplot
! read(2) NZplot
! read(2) Rplot
! read(2) Zplot
! read(2) Zshift
! close(2)
!
! allocate(RePsi(NZplot,NRplot),ImPsi(NZplot,NRplot))
! allocate(Basis(NZplot,NRplot))
!
! open(1, file = 'Solution/NtimesChar',form='unformatted')
! write(*,*) int((to-ti)/dt)
! write(1) 1+int((to-ti)/dt)
! close(1)
!
! open(1, file = 'Solution/timesChar',form='unformatted')
! write(*,*) int((to-ti)/dt)
! write(1) 1+int((to-ti)/dt)
! Do t = ti, to, dt
!
!    call leadingzerostr(t, timeCh)
!    write(*,*) "t",t,timeCh
!    write(1) timeCh
!
!    Basis = c0
!    RePsi = 0.d0
!    ImPsi = 0.d0
!
!    Do L = 0, Lmax
!       write (filenumber2,fmt2) L
!       Do n = 1, Nmax
!          write (filenumber1,fmt1) n
!          Basis = 0.d0
!          open(2, file = 'Basis2D/CoulombBox_L'//trim(filenumber2)//'_n'//trim(filenumber1),form='unformatted')
!          read(2) Basis
!          close(2)
!          RePsi = RePsi + dreal(Basis*Cn_f(n,L,0)*exp(-ci*t*Energias(n,L)))
!          ImPsi = ImPsi + dimag(Basis*Cn_f(n,L,0)*exp(-ci*t*Energias(n,L)))
!
!       enddo
!       Do n = L + 1, Nbounds
!          write (filenumber1,fmt1) n
!          Basis = 0.d0
!          open(2, file = 'Basis2D/BoundBox_L'//trim(filenumber2)//'_n'//trim(filenumber1),form='unformatted')
!          read(2) Basis
!          close(2)
!          RePsi = RePsi + dreal(Basis*Cn_b(n,L,0)*exp(-ci*t*Ebounds(n,L)))
!          ImPsi = ImPsi + dimag(Basis*Cn_b(n,L,0)*exp(-ci*t*Ebounds(n,L)))
!       enddo
!    enddo
!
!    open(2, file = 'Solution/ReSolution_teq_'//timeCh,form='unformatted')
!    write(2) RePsi
!    close(2)
!    open(2, file = 'Solution/ImSolution_teq_'//timeCh,form='unformatted')
!    write(2) ImPsi
!    close(2)
!
! enddo
! close(1)
!
! deallocate(Basis)
! deallocate(RePsi,ImPsi)
!
!
! endsubroutine ComputeFunction
!
! subroutine ExportFunctionGrid
! use ModuleConstants
! use ModuleParameters
! implicit none
! logical             :: dir_e
! real*8, allocatable :: rhopl(:), zpl(:)
! inquire(file='./Solution/test', exist=dir_e)
! if (dir_e) then
! !   write(*,*) "dir exists!",dir_e
! else
!   call system('mkdir -p ./Solution')
!   Open(unit=10,file="Solution/test")
!   close(10)
! end if
!
!
! open(2, file = 'Basis2D/constants',form='unformatted')
! read(2) Nplot
! read(2) NRplot
! read(2) NZplot
! read(2) Rplot
! read(2) Zplot
! read(2) Zshift
! close(2)
!
! allocate(rhopl(NRplot),zpl(NZplot))
! rhopl   = 0.d0
! zpl     = 0.d0
! open(2, file = 'Basis2D/rho',form='unformatted')
! read(2) rhopl
! close(2)
! open(2, file = 'Basis2D/z',form='unformatted')
! read(2) zpl
! close(2)
!
!
! open(2, file = 'Solution/RhoSize.dat')
! write(2,*) NRplot
! close(2)
! open(2, file = 'Solution/ZSize.dat')
! write(2,*) NZplot
! close(2)
! open(2, file = 'Solution/Zshift.dat')
! write(2,*) Zshift
! close(2)
!
!
! open(2, file = 'Solution/RhoGrid',form='unformatted')!, access="stream"
! write(2) rhopl
! close(2)
!
! open(2, file = 'Solution/ZGrid',form='unformatted')
! write(2) zpl
! close(2)
!
! deallocate(rhopl,zpl)
!
! endsubroutine ExportFunctionGrid
!
! subroutine leadingzerostr(x, str_temp)
!     real*8 x
!     character(*) str_temp
!     character(len=8) str
!     character, dimension(3):: signchr = (/'-', ' ', '+' /)
!     write(str,'(F8.4)') 100.0 + abs(x)
!     str_temp = str
!     str_temp(1:1) = signchr(int(sign(1.0,x)) + 2)
! end subroutine leadingzerostr
!


! subroutine ReadBasis2D
! use ModuleConstants
! use ModuleParameters
! implicit none
! integer             :: L, n, ic, id
! Real*8              :: rho2, z2, r
! character(len=8)    :: fmt1, fmt2 ! format descriptor
! character(5)        :: filenumber1 , filenumber2
! real*8,  external   :: BoundCoulomb, FreeCoulomb
! logical             :: dir_e
! complex*16, allocatable :: RePsi(:,:)
! real*8, allocatable :: rhopl(:), zpl(:)
! inquire(file='./Basis2D/test', exist=dir_e)
! if (dir_e) then
! !   write(*,*) "dir exists!",dir_e
! else
!   call system('mkdir -p ./Basis2D')
!   Open(unit=10,file="Basis2D/test")
!   close(10)
! end if
!
!
! open(2, file = 'Basis2D/constants',form='unformatted')
! read(2) Nplot
! read(2) NRplot
! read(2) NZplot
! read(2) Rplot
! read(2) Zplot
! close(2)
!
! allocate(RePsi(NRplot,NZplot),rhopl(NRplot),zpl(NZplot))
! rhopl   = 0.d0
! zpl     = 0.d0
! open(2, file = 'Basis2D/rho',form='unformatted')
! read(2) rhopl
! close(2)
! open(2, file = 'Basis2D/z',form='unformatted')
! read(2) zpl
! close(2)
!
!
! fmt1 = '(I3.3)'
! fmt2 = '(I2.2)'
!
!
! Do L = 0, Lmax
!    write (filenumber2,fmt2) L
!    Do n = 1, Nmax
!       write (filenumber1,fmt1) n
!       RePsi = c0
!       open(2, file = 'Basis2D/CoulombBox_L'//trim(filenumber2)//'_n'//trim(filenumber1),form='unformatted')
!       read(2) RePsi
!       close(2)
!    enddo
!
!    Do n = L+1,Nbounds
!       write (filenumber1,fmt1) n
!       RePsi = c0
!       open(2, file = 'Basis2D/BoundBox_L'//trim(filenumber2)//'_n'//trim(filenumber1),form='unformatted')
!       read(2) RePsi
!       close(2)
!    enddo
! enddo
!
! deallocate(RePsi,rhopl,zpl)
!
!
! endsubroutine ReadBasis2D



