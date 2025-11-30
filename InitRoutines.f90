
module subrutinas

contains

   subroutine angular_archivo()  !!este hace un archivo para guardar los valores de la parte angular indep de Rmax una sola vez
      use funciones_angulares
      use ModuleParameters

      implicit none
      integer :: imax_lpmp, imax_lppmpp, idx_total, idx_lm, idx_lpmp, idx_lppmpp, l, m, lp, mp, lpp, mpp
      real*8 :: res1, res2
      real*8, allocatable :: factor1(:), factor2(:), factor3(:), factor4(:)
      integer :: N_angTE, N_angTM


      integer, allocatable :: block_sizeTE(:), offset_l_TE(:), block_sizeTM(:), offset_l_TM(:)
  
      
      ! N_angTE = ((Lmax+1)**2) * ((Lmax+1)**2) * Lmodes_TE* (Lmodes_TE+2)
      ! N_angTM = ((Lmax+1)**2) * ((Lmax+1)**2) * Lmodes_TM* (Lmodes_TM+2)

      N_angTE = ((Lmax+1)**2) * ((Lmax+1)**2) * ((2*m0+1)*Lmodes_TE - m0*(m0-1))
      N_angTM = ((Lmax+1)**2) * ((Lmax+1)**2) * ((2*m0+1)*Lmodes_TM - m0*(m0-1))
      
      allocate(factor1(N_angTE), factor2(N_angTE), factor3(N_angTM), factor4(N_angTM))

      imax_lpmp   = (Lmax+1)*(Lmax+1)
      imax_lppmpp = (Lmax+1)*(Lmax+1)

      !=== Tamaño de bloque para cada l ===
      allocate(block_sizeTE(Lmodes_TE), block_sizeTM(Lmodes_TM))
      do l = 1, Lmodes_TE
         block_sizeTE(l) = 2 * min(l, m0) + 1
      end do

      do l = 1, Lmodes_TM
         block_sizeTM(l) = 2 * min(l, m0) + 1
      end do

      ! === Desplazamientos offset_l ===
      allocate(offset_l_TE(0:Lmodes_TE), offset_l_TM(0:Lmodes_TM))

      offset_l_TE =0
      offset_l_TE(0) = 0
      do l = 2, Lmodes_TE
         offset_l_TE(l) = offset_l_TE(l-1) + block_sizeTE(l-1)
      end do

      offset_l_TM(1) = 0   ! el bloque l=1 empieza en 0
      offset_l_TM =0
      offset_l_TM(0) = 0
      do l = 2, Lmodes_TM
         offset_l_TM(l) = offset_l_TM(l-1) + block_sizeTM(l-1)
      end do
      offset_l_TM(1) = 0   ! el bloque l=1 empieza en 0


      
      factor1 = 0.d0
      factor2 = 0.d0
      factor3=0.d0
      factor4=0.d0

      !! TE
      do l = 1,Lmodes_TE
         do m = -min(l, m0), min(l, m0)
            idx_lm = offset_l_TE(l) + (m + min(l, m0)) + 1 
            do lp =0, Lmax
               do mp = -lp, lp
                  idx_lpmp = lp*lp + (mp + lp + 1)
                  do lpp = 0, Lmax
                    do mpp =-lpp, lpp
                     idx_lppmpp = lpp*lpp + (mpp + lpp + 1)

                     idx_total = ((idx_lm - 1)*imax_lpmp + (idx_lpmp - 1))*imax_lppmpp + idx_lppmpp
                     call factor_angular(lpp,mpp,l,m,lp,mp,res1) !acompana a Ri
                     call factor_angular(lpp,mpp,l,-m,lp,mp,res2) !para el conjugado
                     factor1(idx_total) = res1
                     factor2(idx_total) = res2*(-1.d0)**(m+1)
                  
                    end do
                  end do
               end do
            end do
         end do
      end do


      !!! TM
      do l = 1,Lmodes_TM
         do m = -min(l, m0), min(l, m0)
            idx_lm = offset_l_TM(l) + (m + min(l, m0)) + 1
            do lp =0, Lmax
               do mp = -lp, lp
                  idx_lpmp = lp*lp + (mp + lp + 1)
                  do lpp = 0, Lmax
                    do mpp =-lpp, lpp
                     idx_lppmpp = lpp*lpp + (mpp + lpp + 1)

                     idx_total = ((idx_lm - 1)*imax_lpmp + (idx_lpmp - 1))*imax_lppmpp + idx_lppmpp
                     
                     factor3(idx_total) = II3(lpp,mpp,l,m,lp,mp) !acompana a R1 y R2
                     factor4(idx_total) = II3(lpp,mpp,l,-m,lp,mp)*(-1.d0)**(m) !para el conjugado

                    end do
                  end do
               end do
            end do
         end do
      end do

      open(unit=34,file="angular/factor_angular1",form='unformatted', status = 'replace')
      write(34) factor1
      close(34)


      open(unit=34,file="angular/factor_angular2",form='unformatted', status = 'replace')
      write(34) factor2
      close(34)

      
      open(unit=34,file="angular/factor_angular3",form='unformatted', status = 'replace')
      write(34) factor3
      close(34)


      open(unit=34,file="angular/factor_angular4",form='unformatted', status = 'replace')
      write(34) factor4
      close(34)

      deallocate(factor1, factor2, factor3, factor4)
      deallocate(block_sizeTE, offset_l_TE, block_sizeTM, offset_l_TM)


   end subroutine angular_archivo


   subroutine leer_angular()
      use ModuleParameters
      implicit none

      logical :: existe
      integer :: N_angTE, N_angTM

      ! N_angTE = ((Lmax+1)**2) * ((Lmax+1)**2) * Lmodes_TE* (Lmodes_TE+2)
      ! N_angTM = ((Lmax+1)**2) * ((Lmax+1)**2) * Lmodes_TM* (Lmodes_TM+2)

      N_angTE = ((Lmax+1)**2) * ((Lmax+1)**2) * ((2*m0+1)*Lmodes_TE - m0*(m0-1))
      N_angTM = ((Lmax+1)**2) * ((Lmax+1)**2) * ((2*m0+1)*Lmodes_TM - m0*(m0-1))

      allocate(factor_angular1(N_angTE), factor_angular2(N_angTE), factor_angular3(N_angTM), factor_angular4(N_angTM))

      inquire(file="angular/factor_angular1", exist=existe)
      if (existe) then
         open(unit=34, file="angular/factor_angular1", form="unformatted", status="old")
         read(34) factor_angular1
         close(34)
      end if
      inquire(file="angular/factor_angular2", exist=existe)
      if (existe) then
         open(unit=34, file="angular/factor_angular2", form="unformatted", status="old")
         read(34) factor_angular2
         close(34)
      end if
      inquire(file="angular/factor_angular3", exist=existe)
      if (existe) then
         open(unit=34, file="angular/factor_angular3", form="unformatted", status="old")
         read(34) factor_angular3
         close(34)
      end if
      inquire(file="angular/factor_angular4", exist=existe)
      if (existe) then
         open(unit=34, file="angular/factor_angular4", form="unformatted", status="old")
         read(34) factor_angular4
         close(34)
      end if

   end subroutine leer_angular

   subroutine leer_matrices()  !!rutina para leer las matrices y dejarlas guardadas en Module parameters
      use ModuleParameters
      implicit none 
      logical :: existe

     ! Si no están asignadas aún, las reservamos
      if (.not. allocated(Txtom2))    allocate(Txtom2(Nfbox, Nfbox))
      if (.not. allocated(Tsurface))  allocate(Tsurface(Nfbox, Nfbox))
      if (.not. allocated(Tcoul))     allocate(Tcoul(Nfbox, Nfbox))


      inquire(file="Hmat/Txtom2", exist=existe)
      if (existe) then
         open(unit=7, file="Hmat/Txtom2", form="unformatted", status="old")
         read(7) Txtom2
         close(7)
      else
         print *, "No se encontró el archivo: Hmat/Txtom2"
      end if

      inquire(file="Hmat/Tsurface", exist=existe)
      if (existe) then
         open(unit=8, file="Hmat/Tsurface", form="unformatted", status="old")
         read(8) Tsurface
         close(8)
      else
         print *, "No se encontró el archivo: Hmat/Tsurface"
      end if


      inquire(file="Hmat/Tcoul", exist=existe)
      if (existe) then
         open(unit=9, file="Hmat/Tcoul", form="unformatted", status="old")
         read(9) Tcoul
         close(9)
      else
         print *, "No se encontró el archivo: Hmat/Tcoul"
      end if


   end subroutine leer_matrices


   subroutine int_auxiliar_material()

      use ModuleParameters
      use pi_module
      use spherical_bessel
      use sph_bessel_roots
      implicit none

      integer :: n, np,i, l, n1, n2, stat, iostat
      real*8  :: func1, func2, func3, func4, prod,x, prod1, prod2
      character(len=2) :: lch, nch
      character(len=8)    :: fmt2
      write(*,*) "int_auxiliar_material"

      fmt2 = '(I2.2)'

      allocate(Txtom2(1:Nfbox,1:Nfbox),Tsurface(1:Nfbox,1:Nfbox),Tcoul(1:Nfbox,1:Nfbox))

      Txtom2    = 0.d0
      Tsurface  = 0.d0
      Tcoul     = 0.d0

      !HAGO LAS INTEGRALES --- elementos de matriz
      allocate(Xquad(1:Nquad), Wquad(1:Nquad))  !son los pesos para el metodo de integracion
      call gauleg(0.d0,1.d0,Xquad,Wquad,Nquad) !los primeros son lo slimites de integracion
      !!estos pesos no dependen de la funcion a integrar, solo de los limites de integracion


      !!independiente de Rmax
      Do n = 1, Nfbox
         Do np = 1, Nfbox
            func1 = 0.d0
            func2 = 0.d0
            func3 = 0.d0
            func4 = 0.d0
            Do i = 1, Nquad
               x = Xquad(i)

               prod = Wquad(i)*sin(np*pi*x)*sin(n*pi*x)*(2.d0)

               func4 = func4 + prod
               func1  = func1   + prod*(1.d0/(x**2.d0)) !termino de barrera centrifuga
               func2  = func2 + prod*(1.d0/(1.d0-(x**2.d0))) !termino de superficie
               func3  = func3    + prod*(1.d0/(x)) !termino de coulomb
            End do

            Txtom2(np,n)    = func1
            Tsurface(np,n)  = func2
            Tcoul(np,n)     = func3

         End Do
      End Do

      open(unit=7, file='Hmat/Txtom2', status='replace', form='unformatted')
      write(7) Txtom2
      close(7)

      open(unit=8,file="Hmat/Tsurface",form='unformatted')
      write(8) Tsurface
      close(8)

      open(unit=9,file="Hmat/Tcoul",form='unformatted')
      write(9) Tcoul
      close(9)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Do l = 1, Lmodes_TE  !barre sobre modos fotonicos
         Do n = 1, Rmodes_TE 
            if (allocated(S)) deallocate(S)
            allocate(S(Nfbox, Nfbox), stat= stat)
            if (stat /= 0) then
               write(*,*) 'Error al hacer allocate S. stat=', stat
               stop
            end if
            S = 0.d0

            Do n1 = 1, Nfbox   !sobre la base de funciones de senos
               Do n2 = 1, Nfbox   

                  func1 = 0.d0
                  Do i = 1, Nquad !!este es para hacer la integral
                     x = Xquad(i)
                     prod1 = Wquad(i)*sph_bessel_j(l,zBessTE(l,n)*x)*sin(n1*pi*x)*sin(n2*pi*x)/x
   
                     func1 = func1 + prod1

                  End Do

                  S(n1,n2) = func1  !!armo las matrices cuyos indices corresponden a los n1 y n2 de la base de funcione senos

               End Do
            End Do


            !!aca ya construyo la matriz para un l y un determinado
            write(lch,fmt2) l
            write(nch, fmt2) n
            open(unit=35, file='Hmat/S_l'//trim(lch)//'_n'//trim(nch), status='replace', form='unformatted', iostat=iostat)
            if (iostat /= 0) then
               write(*,*) 'Error abriendo archivo para l=', l, ' n=', n, ' iostat=', iostat
               stop
            end if
            write(35) S
            close(35)


         
         End Do
      End Do   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Do l = 1, Lmodes_TM  !barre sobre modos fotonicos
         Do n = 1, Rmodes_TM 
            if (allocated(S1)) deallocate(S1)
            allocate(S1(Nfbox, Nfbox))
            if (allocated(S2)) deallocate(S2)
            allocate(S2(Nfbox, Nfbox))

            S1 = 0.d0
            S2 = 0.d0

            Do n1 = 1, Nfbox   !sobre la base de funciones de senos
               Do n2 = 1, Nfbox   

                  func1 = 0.d0
                  func2= 0.d0
                  Do i = 1, Nquad !!este es para hacer la integral
                     x = Xquad(i)
                     prod1 = Wquad(i)*sph_bessel_j(l,zBessTM(l,n)*x)*( -sin(n1*pi*x)*sin(n2*pi*x)/(x**2) + sin(n1*pi*x)*cos(n2*pi*x)*(1.d0/x)*n2*pi ) 
                     ! prod2 = Wquad(i)*sin(n1*pi*x)*sin(n2*pi*x)*(1.d0/x**2)*( sph_bessel_j(l,zBessTM(l,n)*x) + x*dj_l(l,zBessTM(l,n)*x) )
                     prod2 = Wquad(i)*sin(n1*pi*x)*sin(n2*pi*x)*(1.d0/x**2)*( sph_bessel_j(l,zBessTM(l,n)*x) + zBessTM(l,n)*x*dj_l(l,x) )
                     func1 = func1 + prod1
                     func2 = func2 + prod2

                  End Do

                  S1(n1,n2) = func1  !!armo las matrices cuyos indices corresponden a los n1 y n2 de la base de funcione senos
                  S2(n1,n2) = func2
               End Do
            End Do


            !!aca ya construyo la matriz para un l y un determinado
            write(lch,fmt2) l
            write(nch, fmt2) n
            open(unit=36, file='Hmat/S1_l'//trim(lch)//'_n'//trim(nch), status='replace', form='unformatted')
            write(36) S1
            close(36)
            open(unit=37, file='Hmat/S2_l'//trim(lch)//'_n'//trim(nch), status='replace', form='unformatted')
            write(37) S2
            close(37)


         End Do
      End Do 



      ! deallocate(Txtom2, Tsurface, Tcoul)
      deallocate(Xquad, Wquad)

   end subroutine int_auxiliar_material


   subroutine InitMaterialFunctions()
      use ModuleParameters
      use pi_module

      implicit none
      integer             :: n, np, l
      character(len=10) :: Lstate, Nstate
      character(len=8)    :: fmt1,fmt2
      real*8 :: x, func1
      

      real*8, allocatable :: WORK(:)
      integer             :: LWORK, INFO

      write(*,*) "InitMaterialFunctions"


      if (allocated(Hmaterial)) deallocate(Hmaterial)
      if (allocated(WORK)) deallocate(WORK)
      if (allocated(Energias)) deallocate(Energias)
      if (allocated(Estados)) deallocate(Estados)
      if (allocated(Xquad)) deallocate(Xquad)
      if (allocated(Wquad)) deallocate(Wquad)


      allocate(Hmaterial(1:Nfbox,1:Nfbox), stat=info)
      if (info /= 0) then
          write(*,*) "ERROR allocating Hmaterial"
          return
      endif

      
      allocate(Estados(1:Nmax,0:Lmax,1:Nfbox), stat=info)
      if (info /= 0) then
          write(*,*) "ERROR allocating Estados"
          return
      endif
      
      allocate(Energias(1:Nmax,0:Lmax), stat=info)
      if (info /= 0) then
          write(*,*) "ERROR allocating Energias"
          return
      endif


      Estados  = 0.d0
      Energias = 0.d0


      Etilde = -(Zn/Rmax)+((Zn**2.d0)/(2.d0*Rmax))


      LWORK = max(1,3*Nfbox-1)
      allocate(WORK(LWORK))


      fmt1 = '(I3.3)'
      fmt2 = '(I2.2)'



      Do l = 0, Lmax
         if (allocated(En)) deallocate(En)
         allocate(En(1:Nfbox), stat=info)
         if (info /= 0) then
            write(*,*) "ERROR allocating En"
            return
         endif


         write(Lstate,fmt2) l
      !    write(*,*) "Momento angular",l
      !va agragando todos los elementos de matriz para un dado l, aca ya hay dependencia en Rmax
         En        = 0.d0
         Hmaterial = 0.d0
         Hmaterial = ((l*(l+1))/(2.d0*Rmax*Rmax))*Txtom2
         Hmaterial = Hmaterial - (1.d0/(2.d0*Rmax))*Tsurface
         Hmaterial = Hmaterial + (Zn*Ze/Rmax)*Tcoul


         Do n = 1, Nfbox
            Hmaterial(n,n) = Hmaterial(n,n) + 0.5d0*((n*pi/Rmax)**2.d0)  !para la parte diagonal
         End Do


         !!DIAGONALIZACION

         call dsyev('V', 'L', Nfbox, Hmaterial, Nfbox, En, WORK, LWORK, INFO)
         write(*,*) "INFO",INFO

         Do n = 1, Nmax  !!este indice corre sobre el autoestado (autovector) ,columna
            Do np = 1, Nfbox  !!este sobre las componentes de cada autovector, fila
               Estados(n,l,np) = Hmaterial(np,n)
            End do
            Energias(n,l)   = En(n) - Etilde !!!aca estan las energias
         End Do

         
         !!esto es para luego graficar las R(r)
         Do n = 1, Nmax                   !!!!ESTAS SON LA u(r), no las R(r) completas
      !       write(*,*) "En",n,En(n),-Zn*Ze/(2.d0*n*n)
            write(Nstate,fmt1) n
            open(2, file = 'EstadosMateriales/Estado_l'//trim(Lstate)//'_n'//trim(Nstate)//'.dat', status = 'replace')

            Do x = Rmax/1000.d0, Rmax, Rmax/5000.d0
               func1 = 0.d0
               Do np = 1, Nfbox
                  func1 = func1 + Hmaterial(np,n)*sin(np*pi*x/Rmax)*((2.d0/Rmax)**0.5d0) !!Hmaterial(np,n) tiene los coeficientes del autovector n y barre sobre np que va hasta Nfbox
               End Do
               Write(2,*) x, func1
            End do
            close(2)
         End Do

      End do  !!--este End Do es del Do que va sobre los valores de l

      deallocate(WORK,Hmaterial,En)



   End subroutine InitMaterialFunctions


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine InitZerosBessel
      use ModuleParameters
      use pi_module
      use sph_bessel_roots
      use spherical_bessel
      implicit none
      integer :: l,n

      write(*,*) "InitBessel"   !!sacar a eso dek loop prob con allocate

      allocate(zBessTE(Lmodes_TE,Rmodes_TE)) !!sus filas son todos los ceros de bessel para ese l dado
      allocate(zBessTM(Lmodes_TM,Rmodes_TM))


      zBessTE = 0.d0
      zBessTM = 0.d0


      Do l=1, Lmodes_TE
         Do n = 1,Rmodes_TE
            zBessTE(l,n) = bessel_j_root(l,n)
         End Do
      End Do

      Do l=1, Lmodes_TM
         Do n = 1,Rmodes_TM
            zBessTM(l,n) = bessel_jprime_root(l,n)
         End Do
      End Do

      !!!!CHECKEADO QUE LAS RAICES ESTAN BIEN
   End subroutine InitZerosBessel



   subroutine ctes_normalizacion()  !!dependen de Rmax

      use ModuleParameters
      use sph_bessel_roots
      use spherical_bessel
      implicit none
      integer :: l,n
      real*8 :: xln, yln


      
      allocate(NTE(Lmodes_TE,Rmodes_TE))
      allocate(NTM(Lmodes_TM,Rmodes_TM))
      NTE = 0.d0
      NTM = 0.d0

      !!!!factor de normalizacion de los modos TE y TM, dependen de Rmax
      Do l=1, Lmodes_TE
         Do n = 1,Rmodes_TE
            xln = zBessTE(l,n)
            NTE(l,n) = 1.d0/sqrt(0.5d0 * (Rmax**3) * (dj_l(l,xln))**2 )
         End Do
      End Do

      Do l=1, Lmodes_TM
         Do n = 1,Rmodes_TM
            yln = zBessTM(l,n)
            NTM(l,n) = 1.d0 / sqrt(0.5d0 * (Rmax**3) * (ddxj(l,yln)**2)/( yln**2 - l*(l+1) ))
         End Do
      End Do

   end subroutine ctes_normalizacion


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! subroutine InitRadialIntegrals
   !    use ModuleParameters
   !    use pi_module
   !    use sph_bessel_roots
   !    use spherical_bessel
   !    use ModuleState
   !    implicit none

   !    integer  :: l, m, n, np, lp, mp, lpp, mpp, npp, n1, n2
   !    integer  :: i, ip, ipp, icampo
   !    real*8   :: x, prod1, prod2, func1, func2, func3, func4
   !    real*8, allocatable :: SinCosBess1(:,:), SinCosBess2(:,:), vec1(:), vec2(:), vec3(:), vec4(:)
   !    character(2)        :: Nch, Lch, Mch, icampoch
   !    character(2)        :: ModphCh
   !    character(len=8)    :: fmt1,fmt2

   !    character(len=200) :: filename_plus, filename_minus, filenameTE, filenameTE_CONJ
   !    character(len=3)   :: llch, nnch
   !    integer :: irow,icol
   !    real*8, allocatable :: RiModeTE(:,:), RiModeTEconj(:,:)


   !    write(*,*) "InitRadialIntegrals"

   !    fmt1 = '(I3.3)'
   !    fmt2 = '(I2.2)'




   !    ! n = 1
   !    ! write(Nch,fmt2) n
   !    ! l = 1
   !    ! write(Lch,fmt2) l
   !    ! open(unit=10,file="matrices/RiPlus_L"//Lch//"_n"//Nch,form='unformatted')
   !    ! read(10) RiPlus
   !    ! close(10)



   !    !HAGO LAS INTEGRALES --- elementos de matriz
   !    allocate(Xquad(1:Nquad), Wquad(1:Nquad))  !son los pesos para el metodo de integracion
   !    call gauleg(0.d0,1.d0,Xquad,Wquad,Nquad) !los primeros son lo slimites de integracion
   !    !!estos pesos no dependen de la funcion a integrar, solo de los limites de integracion


   !    allocate(RiPlus(NsysMA,NsysMA), RiMinus(NsysMA,NsysMA)) !--> son las matrices que van a ir en los bloques para un dado modo fotonico. Esa matriz ira en la matriz grande si es que las deltas lo permiten
   !    allocate(SinCosBess1(Nfbox,Nfbox),SinCosBess2(Nfbox,Nfbox))
   !    allocate(vec1(Nfbox), vec2(Nfbox), vec3(Nfbox), vec4(Nfbox))
   !    RiPlus = 0.d0
   !    RiMinus = 0.d0

   !    allocate(RiModeTE(NsysMA,NsysMA), RiModeTEconj(NsysMA,NsysMA))
   !    RiModeTE = 0.d0
   !    RiModeTEconj = 0.d0
   !    ! Do en angulares de campo y materiales para


   !    Do l = 1, Lmodes_TE  !barre sobre modos fotonicos
   !       Do n = 1, Rmodes_TE 

   !          RiPlus  = 0.d0
   !          RiMinus = 0.d0


   !          SinCosBess1 = 0.d0
   !          SinCosBess2 = 0.d0
   !          Do n1 = 1, Nfbox   !sobre la base de funciones de senos
   !             Do n2 = 1, Nfbox   

   !                func1 = 0.d0
   !                func2 = 0.d0
   !                Do i = 1, Nquad !!este es para hacer la integral
   !                   x = Xquad(i)
   !                   prod1 = Wquad(i)*(2.d0)*sph_bessel_j(l,zBessTE(l,n)*x)*sin(n1*pi*x)*cos(n2*pi*x)*(n2*pi/Rmax)
   !                   prod2 = Wquad(i)*(2.d0)*sph_bessel_j(l,zBessTE(l,n)*x)*sin(n1*pi*x)*(1.d0/(Rmax*x))*sin(n2*pi*x)
   !                   func1 = func1 + prod1
   !                   func2 = func2 + prod2
   !                End Do
   !                SinCosBess1(n1,n2) = func1  !!armo las matrices cuyos indices corresponden a los n1 y n2 de la base de funcione senos
   !                SinCosBess2(n1,n2) = func2
   !             End Do
   !          End Do


   !          Do lpp = 0, Lmax
   !             Do npp = 1, Nmax

   !                vec1(:)= Estados(npp,lpp,:) !estos son los a_np que expanden en la base de senos a cada autovector nl
   !    !             ip = 0
   !                Do lp = 0, Lmax
   !                   Do np = 1, Nmax
   !    !                   ip = ip + 1

   !                      vec2(:) = Estados(np,lp,:)

   !                      vec3(:) = matmul(SinCosBess1,vec2)

   !                      func1 = DOT_PRODUCT(vec1,vec3)

   !                      vec4(:) = matmul(SinCosBess2,vec2)

   !                      func2 = DOT_PRODUCT(vec1,vec4)

   !                      prod1 = func1 +lp*func2  !Riplus
   !                      prod2 = func1 -(lp+1)*func2 !!Riminus

   !                      Do mpp = -lpp, lpp
   !                         ipp = (lpp*(lpp+1) + mpp)*Nmax + npp
   !                         Do mp = -lp, lp
   !                            ip = (lp*(lp+1) + mp)*Nmax + np
   !                            RiPlus(ipp,ip)  = prod1                 !!aca estoy llenando las matrices con los valores repetidos para mp y mpp, pues no dependen de ellos
   !                            RiMinus(ipp,ip) = prod2
   !                            write(*,'(A,I0,A,I0,A,3I0,A,3I0)') '(ipp,ip)=', ipp, ',', ip, ' lpp,mpp,npp=', lpp, mpp, npp, ' lp,mp,np=', lp, mp, np
   !                            !estan bien asignados los indices
                           

   !                         End do
   !                      End Do

   !                   End do
   !                End Do
   !             End do
   !          End do


   !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !    !!para un checkeo 



   !          ! Convertir l y n a cadenas con formato fijo (por ejemplo l02 n03)
   !          write(llch, '(I2.2)') l
   !          write(nnch, '(I3.2)') n

   !          ! Construir los nombres de archivo
   !          filename_plus  = "matrices/RiPlus_l"  // trim(llch) // "_n" // trim(nnch) // ".dat"
   !          filename_minus = "matrices/RiMinus_l" // trim(llch) // "_n" // trim(nnch) // ".dat"

   !          ! Abrir y escribir en formato texto (replace: sobreescribe si ya existe)
   !          open(unit=10, file=filename_plus,  status='replace', action='write', form='formatted')
   !          open(unit=11, file=filename_minus, status='replace', action='write', form='formatted')

   !          ! Guardar matrices (usa bucles si querés formato legible)
   !          do irow = 1, size(RiPlus,1)
   !             write(10, '(1000ES20.10)') (RiPlus(irow,icol), icol=1,size(RiPlus,2))
   !          end do
   !          do irow = 1, size(RiMinus,1)
   !             write(11, '(1000ES20.10)') (RiMinus(irow,icol), icol=1,size(RiMinus,2))
   !          end do

   !          close(10)
   !          close(11)


   !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !          !INTEGRALES SOLO RADIALES CHECKEADAS, DAN CORRECTO      
   !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !          Do m = -l, l


   !             Do lp = 0, Lmax
   !                Do mp = -lp, lp
   !                   Do lpp = 0, Lmax
   !                      Do mpp = -lpp, lpp

   !                         call evaluateangular1(lpp,-mpp,l,m,lp,mp,func1) !acompana a Riminus
   !                         call evaluateangular2(lpp,-mpp,l,m,lp,mp,func2) !  acompana a Riplus

   !                         call evaluateangular1(lpp,-mpp,l,-m,lp,mp,func3) !
   !                         call evaluateangular2(lpp,-mpp,l,-m,lp,mp,func4) !
   !                         !!!el -mpp es porque la funcion '' va conjugada --> el armonico esfercio se multiplica por un (-1)^mpp y mpp->-mpp
   !                         func1 = func1*(-1.d0)**(mpp)
   !                         func2 = func2*(-1.d0)**(mpp)
   !                         func3 = func3*(-1.d0)**(mpp)
   !                         func4 = func4*(-1.d0)**(mpp)

   !                         Do np = 1, Nmax
   !                            ip = (lp*(lp+1) + mp)*Nmax + np
   !                            Do npp = 1, Nmax
   !                               ipp = (lpp*(lpp+1) + mpp)*Nmax + npp

   !                               RiModeTE(ipp,ip)=RiMinus(ipp,ip)*func1+RiPlus(ipp,ip)*func2

   !                               RiModeTEconj(ipp,ip)=(RiMinus(ipp,ip)*func3+RiPlus(ipp,ip)*func4)*((-1.d0)**(m+1))


   !                            End Do
   !                         End Do
   !                      End Do
   !                   End Do
   !                End Do
   !             End Do


   !             icampo = (l*l -1 + m + l)*Rmodes_TE + n  !!este indice esta bien si es desde l=1, no l=0
   !             write(icampoch,fmt2) icampo
   !             open(unit=10,file="matrices/RiPlusMinusTE_Icampo"//icampoch,form='unformatted')
   !             write(10) RiModeTE
   !             close(10)
   !             open(unit=10,file="matrices/RiPlusMinusConjTE_Icampo"//icampoch,form='unformatted')
   !             write(10) RiModeTEconj
   !             close(10)
   !          End Do

   !       End do
   !    End do

   !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! checkeo de parte angular

   !       ! do m = 1, M_TE
   !       !    write(icampoch,fmt2) m
   !       !    filenameTE  = "matrices/RiPlusMinusTE_Icampo" // trim(icampoch) // ".dat"
   !       !    filenameTE_CONJ = filenameTE = "matrices/RiPlusMinusConjTE_Icampo" // trim(icampoch) // ".dat"

   !       !    ! Abrir y escribir en formato texto (replace: sobreescribe si ya existe)
   !       !    open(unit=12, file=filenameTE,  status='replace', action='write', form='formatted')
   !       !    open(unit=13, file=filenameTE_CONJ, status='replace', action='write', form='formatted')

   !       !    ! Guardar matrices (usa bucles si querés formato legible)
   !       !    do irow = 1, size(RiModeTE,1)
   !       !       write(12, '(1000ES20.10)') (RiModeTE(irow,icol), icol=1,size(RiModeTE,2))
   !       !    end do
   !       !    do irow = 1, size(RiModeTEconj,1)
   !       !       write(13, '(1000ES20.10)') (RiModeTEconj(irow,icol), icol=1,size(RiModeTEconj,2))
   !       !    end do

   !       !    close(12)
   !       !    close(13)

   !       ! end do

   !          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   !    deallocate(RiPlus, RiMinus,SinCosBess1,SinCosBess2,vec1,vec2,vec3)
   !    deallocate(Xquad, Wquad)

   !    deallocate(RiModeTE,RiModeTEconj)

   ! End subroutine InitRadialIntegrals



   ! subroutine BuildTotalMatrix()
   !    use complex_module
   !    use ModuleParameters
   !    use pi_module
   !    use sph_bessel_roots
   !    use ModuleState
   !    use clight_module
   !    implicit none
   !    integer :: i,j , l, n, m, imi, ima, jmi, jma
   !    real*8 :: b_a, b_adaga, wln
   !    integer, allocatable :: vector_fock(:), estado_mat(:), mm(:)
   !    integer, allocatable :: fock_izq(:), mat_izq(:),fock_der(:), mat_der(:)


   !    complex*16, allocatable :: Hint(:,:)
   !    real*8, allocatable :: Hm(:), Hf(:) !--> Hamiltonianos
   !    complex*16, allocatable :: H(:,:)  !!--> matriz luego de diagonalizar, sus columnas son los autovectores
   !    real*8, allocatable :: E(:) !!--> vector de autoenergias

   !    real*8 :: E0 !energia de vacio de fotones
   !    character(2)        :: icampoch
   !    character(len=8)    :: fmt1,fmt2

   !    real*8, allocatable :: RiModeTE(:,:), RiModeTEconj(:,:)

   !    integer ::LDA, LWORK, INFO
   !    complex*16, allocatable :: VR(:,:), VL(:,:) !autovectoresderechos
   !    complex*16, allocatable :: W(:)     !autovalores en principio complejos
   !    complex*16, allocatable :: WORK(:) 


   !    fmt1 = '(I3.3)'
   !    fmt2 = '(I2.2)'


   !    allocate(Hint(Nsystem,Nsystem), Hf(Nsystem), Hm(Nsystem))
   !    Hint = c0
   !    Hm = 0.d0
   !    Hf = 0.d0

   !    write(*,*) "BuildTotalMatrix"

   !    ! call InitConfSpace()

   !    !H de fotones
   !    !Energia de una configuracion de fotones
   !    E0 = 0.d0

   !    Do i = 1, M_TE
   !       l = modes_matrix(i,2) !es el l en la matriz de estados fotonicos
   !       n = modes_matrix(i,4) !es el n en la matriz de estados fotonicos
   !       E0 = E0 + bessel_j_root(l, n)
   !    End do


   !    !falta definir dj_j_root para que encuentre las raices para los TM
   !    ! Do i = M_TE + 1, M_TM
   !    !    l = modes_matrix(i,2) !es el l en la matriz de estados fotonicos
   !    !    n = modes_matrix(i,4) !es el n en la matriz de estados fotonicos
   !    !    E0 = E0 + dj_j_root(l, n)
   !    ! End do

   !    E0 = E0*(1.d0/Rmax)*clight*0.5d0

   !    allocate(vector_fock(M_TE+M_TM),estado_mat(3))

   !    !!PARTE DIAGONAL
   !    Do i = 1, Nsystem

   !       call get_state(i,NsysMA,vector_fock,estado_mat) !aca ya tengo los estados fotonico y materiales

   !       Hf(i) = E_N(vector_fock) + E0 !construyo la matriz diagonal de las energias de los fotones (es un vector directamente)

   !       !lmn materiales
   !       l = estado_mat(1)
   !       n = estado_mat(3)
   !       Hm(i) = Energias(n,l) !!cada estado i-esimo viene con su Emat dada por get_state

   !       Hint(i,i) = Hint(i,i) + c1*Hm(i) + c1*Hf(i) !!leno la parte diagonal en Hint

   !    End do



   !    allocate(fock_izq(M_TE + M_TM), mat_izq(3))
   !    allocate(fock_der(M_TE + M_TM), mat_der(3))
   !    allocate(mm(4)) !--> aca estan los conjuntos de indices para cada estado fotonico

   !    allocate(RiModeTE(NsysMA,NsysMA),RiModeTEconj(NsysMA,NsysMA))

   !    !H DE INTERACCION
   !    Do i = 1, NsysPH
   !       imi = (i-1)*NsysMA+1 !este es el inicio de las filas del bloque para esa conf de fock
   !       ima = i*NsysMA       !este es el final de las filas para ese bloque
   !       call get_state(imi, NsysMA,fock_izq, mat_izq )
   !       ! write(*,*)'fila',imi
   !       ! write(*,*)'estado fotonico de fock de izquerda(fila)', fock_izq  

   !       Do j = 1, NsysPH
   !          jmi = (j-1)*NsysMA+1 !este es el inicio de las columnas del bloque para esa conf de fock
   !          jma = j*NsysMA       !este es el inicio de las columanas para ese bloque
   !          call get_state(jmi, NsysMA,fock_der, mat_der)
   !          ! write(*,*)'columna',jmi
   !          ! write(*,*)'estado fotonico de fock de derecha(columna)', fock_der 

   !          Do m = 1,M_TE!+M_TM  !!este barre sobre modos fotonicos dentro de un dado bloque, en vez de checkear todos los modos, veo hasta el primero que no es nulo
   !             !checkeo si son cero o no para un dado m
   !             !obtengo el l correspondiente
   !             l = modes_matrix(m,2)
   !             n = modes_matrix(m,4)
   !             b_a = delta_a(fock_der, fock_izq, m) !destruccion
   !             b_adaga = delta_adaga(fock_der, fock_izq, m) !creacion

   !             ! if (b_a /= 0 .or. b_adaga /= 0) then
   !             !    write(*,*)'f_izq, f_der, m ', fock_izq, '//',fock_der,'//', m
   !             ! end if

   !                write(icampoch,fmt2) m
   !                ! write(*,*) 'ESTADO m: ',modes_matrix(m,:)

   !                open(unit=10,file="matrices/RiPlusMinusTE_Icampo"//icampoch,form='unformatted')
   !                read(10) RiModeTE
   !                close(10)
   !                open(unit=10,file="matrices/RiPlusMinusConjTE_Icampo"//icampoch,form='unformatted')
   !                read(10) RiModeTEconj
   !                close(10)

   !                !frecuencias
   !                wln = bessel_j_root(l,n)*(1.d0/Rmax)*clight

   !                Hint(imi:ima,jmi:jma) = Hint(imi:ima,jmi:jma) + sqrt(1.d0/(2*wln))*sqrt(4*pi/3.d0)*(1.d0/sqrt(l*(l+1)*1.d0))*( -ci*(1.d0/(clight*me))*RiModeTE*b_a - ci*(1/clight*me)*RiModeTEconj*b_adaga) !!aca pega esos bloques en la matriz grande en los rangos indicados
   !                ! exit !!esto es porque al primero distinto de cero, ya se que el resto sera cero

   !          End do

   !       End do
   !    End do
   !    deallocate(RiModeTE,RiModeTEconj)

   !    !!COMENTARIO
   !    !!para una dada configuracion de N y N' hay un solo entre a destruccion y adaga creacion que no va a ser cero, el otro es cero
   !    !! y ademas para solo un modo fotonico m es el que va a contribuir para esa dada conf de N y N'. Si uno es distinto de cero es porque en ese lugar m 
   !    !!se cumplen las condiones, en los otros automaticamente no se cumplen
   !    Do i = 1, Nsystem
   !       Do j= i, Nsystem
   !       write(*,*) "*******************************************"
   !          write(*, *) i,j
   !          write(*, *) Hint(i,j)
   !          write(*, *) Hint(j,i)
   !    !         write(*,  '(2I5, 2F8.3, F8.3, 2F8.3, F8.3)') i,j,dreal(Hint(i,j)),dimag(Hint(i,j)),dreal(Hint(j,i)),dimag(Hint(j,i))
   !       End Do
   !    End Do




   !    stop

   !    !!!DIAG sin asumir que es hermitiana, para checkear
   !    LDA = Nsystem

   !    allocate(W(Nsystem), VR(Nsystem,Nsystem), VL(Nsystem,Nsystem))

   !    LWORK = -1
   !    allocate(WORK(1))
   !    call zgeev('N','V', Nsystem, Hint, LDA, W, VL, LDA, VR, LDA, WORK, LWORK, INFO)

   !    LWORK = int(real(WORK(1)))
   !    deallocate(WORK)
   !    allocate(WORK(LWORK))

   !    !DIAGINALIZACION
   !    call zgeev('N','V', Nsystem, Hint, LDA, W, VL, LDA, VR, LDA, WORK, LWORK, INFO)

   !    if (INFO == 0) then
   !       print *, 'Autovalores:'
   !       print *, W
   !       print *, 'Autovectores derechos (columnas de VR):' !!los autovectores derechos son los usuales, los izq son para la dagada
   !       print *, VR
   !    else
   !       print *, 'Error en diagonalización, INFO=', INFO
   !    end if

   !    ! !H DE INTERACCION
   !    ! Hint = 0.d0
   !    ! Do i = 1, Nsystem
   !    !    get_state(i, NsysMA,fock_izq, mat_izq )
   !    !    Do j = 1, Nsystem
   !    !       get_state(j, NsysMA,fock_der, mat_der)
   !    !       Do m = 1,M_TE+M_TM  !!este suma sobre modos fotonicos
   !    !
   !    !          !checkeo si son cero o no para un dado m
   !    !          b = delta(fock_izq, fock_der, m)
   !    !          if ( b==0 ) then
   !    !             Hint(i,j) = 0
   !    !          else
   !    !
   !    !
   !    !             !calculo el snadwich para ese m = alpha,l,m,n fotonicos y mat_izq,mat_der
   !    !             mm = modes_matrix(m,:)
   !    !             K  = funcion_sandwich(mm, mat_der, mat_izq)
   !    !             Hint(i,j) = Hint(i,j) + K*b  !!--> es el termino para un dado m
   !    !
   !    !
   !    !          end if
   !    !
   !    !       End do
   !    !    End do
   !    ! End do


   ! End subroutine BuildTotalMatrix



   !!!!!!!!!!!! version 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine integral_radial()
      use ModuleParameters
      use pi_module
      use sph_bessel_roots
      use spherical_bessel
      use ModuleState
      use funciones_angulares
      implicit none

      integer  :: l, m, n, np, lp, mp, lpp, mpp, npp, imax_lpmp, imax_lppmpp, idx_total, idx_lm, idx_lpmp, idx_lppmpp
      integer  :: ip, ipp, icampo, offset
      real*8   :: func, res1, res2, func1, func2, res3, res4

      character(2)        :: icampoch, lch, nch
  

      character(len=8)    :: fmt1,fmt2

      real*8, allocatable :: RiaTE(:,:), RiaTE_conj(:,:), RiaTM(:,:), RiaTM_conj(:,:) ,Ri(:,:), R1(:,:), R2(:,:)
      real*8, allocatable :: vec1(:), vec2(:), vec3(:), vec4(:)

      integer, allocatable :: block_sizeTE(:), offset_l_TE(:), block_sizeTM(:), offset_l_TM(:)


      write(*,*) "InitRadialIntegrals"

      fmt1 = '(I3.3)'
      fmt2 = '(I2.2)'

      !!!!!!!!!!!!!!!! PARTE ELECTRICA

      ! !HAGO LAS INTEGRALES --- elementos de matriz
      ! allocate(Xquad(1:Nquad), Wquad(1:Nquad))  !son los pesos para el metodo de integracion
      ! call gauleg(0.d0,1.d0,Xquad,Wquad,Nquad) !los primeros son lo slimites de integracion
      ! !!estos pesos no dependen de la funcion a integrar, solo de los limites de integracion

      if(allocated(Ri)) deallocate(Ri)
      allocate(Ri(NsysMA,NsysMA)) !--> son las matrices que van a ir en los bloques para un dado modo fotonico. Esa matriz ira en la matriz grande si es que las deltas lo permiten
      
      if(allocated(R1) .and. allocated(R2)) deallocate(R1,R2)
      allocate(R1(NsysMA,NsysMA),R2(NsysMA,NsysMA))


      if(allocated(vec1)) deallocate(vec1)
      if(allocated(vec2)) deallocate(vec2)
      if(allocated(vec3)) deallocate(vec3)
      if(allocated(vec4)) deallocate(vec4)
      allocate(vec1(Nfbox), vec2(Nfbox), vec3(Nfbox), vec4(Nfbox))

      if(allocated(RiaTE)) deallocate(RiaTE)
      if(allocated(RiaTE_conj)) deallocate(RiaTE_conj)

      if(allocated(RiaTM)) deallocate(RiaTM)
      if(allocated(RiaTM_conj)) deallocate(RiaTM_conj)


      allocate(RiaTE(NsysMA,NsysMA), RiaTE_conj(NsysMA,NsysMA))
      allocate(RiaTM(NsysMA,NsysMA), RiaTM_conj(NsysMA,NsysMA))


      !! para poder armar los indices para leer la parte angular!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(block_sizeTE(Lmodes_TE), block_sizeTM(Lmodes_TM))
      do l = 1, Lmodes_TE
         block_sizeTE(l) = 2 * min(l, m0) + 1
      end do
     
      do l = 1, Lmodes_TM
         block_sizeTM(l) = 2 * min(l, m0) + 1
      end do

      ! === Desplazamientos offset_l ===
      allocate(offset_l_TE(0:Lmodes_TE), offset_l_TM(0:Lmodes_TM))

      offset_l_TE =0
      offset_l_TE(0) = 0
      do l = 2, Lmodes_TE
         offset_l_TE(l) = offset_l_TE(l-1) + block_sizeTE(l-1)
      end do

      offset_l_TM(1) = 0   ! el bloque l=1 empieza en 0
      offset_l_TM =0
      offset_l_TM(0) = 0
      do l = 2, Lmodes_TM
         offset_l_TM(l) = offset_l_TM(l-1) + block_sizeTM(l-1)
      end do
      offset_l_TM(1) = 0   ! el bloque l=1 empieza en 0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      RiaTE = 0.d0
      RiaTE_conj = 0.d0
      RiaTM = 0.d0
      RiaTM_conj = 0.d0

      func= 0.d0
      func1= 0.d0
      func2= 0.d0

      imax_lpmp   = (Lmax+1)*(Lmax+1)
      imax_lppmpp = (Lmax+1)*(Lmax+1)

      !!TE
      icampo =0
      Do l = 1, Lmodes_TE  !barre sobre modos fotonicos
         Do n = 1, Rmodes_TE 

            Ri = 0.d0
            !!solo leo las matrices pero no hago las integrales
            if(allocated(S)) deallocate(S)
            allocate(S(Nfbox,Nfbox))
            write(lch,fmt2) l
            write(nch, fmt2) n
            open(unit=35, file='Hmat/S_l'//trim(lch)//'_n'//trim(nch), form='unformatted')
            read(35) S
            close(35)

            S = S*(2/Rmax)  !!aca corrigo para el Rmax particular



            !!!ver la forma de que esta rutina tome los estados que se van a ir generando a cada R, ver si lo hace como argumento o 
            Do lpp = 0, Lmax
               Do npp = 1, Nmax
                  !!Estados como variable alocatable se guarda en parameters y se usa en cada iteracion en todas las rutinas antes de actualizarse
                  vec1(:)= Estados(npp,lpp,:) !estos son los a_np que expanden en la base de senos a cada autovector nl

                  Do lp = 0, Lmax
                     Do np = 1, Nmax


                        vec2(:) = Estados(np,lp,:)

                        vec3(:) = matmul(S,vec2)

                        func = DOT_PRODUCT(vec1,vec3)  !Ri


                        Do mpp = -lpp, lpp
                           ipp = (lpp*(lpp+1) + mpp)*Nmax + npp
                           Do mp = -lp, lp
                              ip = (lp*(lp+1) + mp)*Nmax + np
                              Ri(ipp,ip)  = func          !!aca estoy llenando las matrices con los valores repetidos para mp y mpp, pues no dependen de ellos
         
                           End do
                        End Do

                     End do
                  End Do
               End do
            End do



            do m = -min(l, m0), min(l, m0)
            idx_lm = offset_l_TE(l) + (m + min(l, m0)) + 1 


               Do lp = 0, Lmax
                  Do mp = -lp, lp
                     idx_lpmp = lp*lp + (mp + lp + 1)
                     Do lpp = 0, Lmax
                        Do mpp = -lpp, lpp
                           idx_lppmpp = lpp*lpp + (mpp + lpp + 1)
                           
                           idx_total = ((idx_lm - 1)*imax_lpmp + (idx_lpmp - 1))*imax_lppmpp + idx_lppmpp  !!Este es un indice global que ordena la tupla (l,m,lp,mp,lpp,mpp)

                           
                           res1 = factor_angular1(idx_total) !!estos son vectores que ya estan leidos una vez al principio y quedan en la memoria
                           res2 = factor_angular2(idx_total)

                           Do np = 1, Nmax
                              ip = (lp*(lp+1) + mp)*Nmax + np
                              Do npp = 1, Nmax
                                 ipp = (lpp*(lpp+1) + mpp)*Nmax + npp

                                 func = Ri(ipp,ip)*res1
                                 if(dabs(func).gt.0.0000000001d0) RiaTE(ipp,ip)      = func
                                 func = Ri(ipp,ip)*res2
                                 if(dabs(func).gt.0.0000000001d0) RiaTE_conj(ipp,ip) = func


                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do

               ! offset = 0
               ! do lp = 1, l-1
               !    offset = offset + (2*lp + 1) * Rmodes_TE
               ! end do
               ! icampo = offset + (m + l) * Rmodes_TE + n
               offset = 0
               do lp = 1, l-1
                  offset = offset + (2*min(lp, m0) + 1) * Rmodes_TE
               end do

               icampo = offset + ( (m + min(l, m0)) * Rmodes_TE ) + n

               ! icampo = (l*l -1 + m + l)*Rmodes_TE + n  !!este indice esta bien si es desde l=1, no l=0
               write(icampoch,fmt2) icampo
               open(unit=15,file="matrices/RiaTE_Icampo"//trim(icampoch),form='unformatted', status = 'replace')
               write(15) RiaTE
               close(15)
               open(unit=15,file="matrices/RiaTEconj_Icampo"//trim(icampoch),form='unformatted', status = 'replace')
               write(15) RiaTE_conj
               close(15)
            End Do

         End do
      End do




   !!!!! PARTE MAGNETICA!!!!!!!!!!!!!!!!!!!!!!!!!

      icampo =0
      Do l = 1, Lmodes_TM  !barre sobre modos fotonicos
         Do n = 1, Rmodes_TM 

            R1 = 0.d0
            R2 = 0.d0
            !!solo leo las matrices pero no hago las integrales
            if(allocated(S1)) deallocate(S1)
            if(allocated(S2)) deallocate(S2)           
            allocate(S1(Nfbox,Nfbox))
            allocate(S2(Nfbox,Nfbox))
            write(lch,fmt2) l
            write(nch, fmt2) n
            open(unit=36, file='Hmat/S1_l'//trim(lch)//'_n'//trim(nch), form='unformatted')
            read(36) S1
            close(36)
            open(unit=37, file='Hmat/S2_l'//trim(lch)//'_n'//trim(nch), form='unformatted')
            read(37) S2
            close(37)

            S1= S1*(2/Rmax**2)  !!aca corrigo para el Rmax particular
            S2= S2*(2/Rmax**2)



            Do lpp = 0, Lmax
               Do npp = 1, Nmax
                  !!Estados como variable alocatable se guarda en parameters y se usa en cada iteracion en todas las rutinas antes de actualizarse
                  vec1(:)= Estados(npp,lpp,:) !estos son los a_np que expanden en la base de senos a cada autovector nl

                  Do lp = 0, Lmax
                     Do np = 1, Nmax


                        vec2(:) = Estados(np,lp,:)

                        vec3(:) = matmul(S1,vec2)
                        func1 = DOT_PRODUCT(vec1,vec3)  !R1

                        vec4(:) = matmul(S2,vec2)       !R2
                        func2 = DOT_PRODUCT(vec1,vec4)


                        Do mpp = -lpp, lpp
                           ipp = (lpp*(lpp+1) + mpp)*Nmax + npp
                           Do mp = -lp, lp
                              ip = (lp*(lp+1) + mp)*Nmax + np
                              R1(ipp,ip)  = func1          !!aca estoy llenando las matrices con los valores repetidos para mp y mpp, pues no dependen de ellos
                              R2(ipp,ip) = func2
                           End do
                        End Do

                     End do
                  End Do
               End do
            End do



            do m = -min(l, m0), min(l, m0)
            idx_lm = offset_l_TM(l) + (m + min(l, m0)) + 1 


               Do lp = 0, Lmax
                  Do mp = -lp, lp
                     idx_lpmp = lp*lp + (mp + lp + 1)
                     Do lpp = 0, Lmax
                        Do mpp = -lpp, lpp
                           idx_lppmpp = lpp*lpp + (mpp + lpp + 1)
                           
                           idx_total = ((idx_lm - 1)*imax_lpmp + (idx_lpmp - 1))*imax_lppmpp + idx_lppmpp !!Este es un indice global que ordena la tupla (l,m,lp,mp,lpp,mpp)

                           
                           res3 = factor_angular3(idx_total) !!estos son vectores que ya estan leidos una vez al principio y quedan en la memoria
                           res4 = factor_angular4(idx_total)

                           Do np = 1, Nmax
                              ip = (lp*(lp+1) + mp)*Nmax + np
                              Do npp = 1, Nmax
                                 ipp = (lpp*(lpp+1) + mpp)*Nmax + npp

                                 RiaTM(ipp,ip)= res3*( (l*(l+1))*R1(ipp,ip) + 0.5d0*( (l*(l+1)) + (lp*(lp+1)) - (lpp*(lpp+1)))*R2(ipp,ip)   )
                                 RiaTM_conj(ipp,ip) = res4*( (l*(l+1))*R1(ipp,ip) + 0.5d0*( (l*(l+1)) + (lp*(lp+1)) - (lpp*(lpp+1)))*R2(ipp,ip)   )


                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do

               
               ! offset = 0
               ! do lp = 1, l-1
               !    offset = offset + (2*lp + 1) * Rmodes_TM
               ! end do
               ! icampo = offset + (m + l) * Rmodes_TM + n

               offset = 0
               do lp = 1, l-1
                  offset = offset + (2*min(lp, m0) + 1) * Rmodes_TM
               end do

               icampo = offset + ( (m + min(l, m0)) * Rmodes_TM ) + n

               ! icampo = (l*l -1 + m + l)*Rmodes_TM + n  !!este indice esta bien si es desde l=1, no l=0
               write(icampoch,fmt2) icampo + M_TE  !!!esto es para que le ponga el nombre despues de los M_TE, ej: si M_TE = 3, que le asigne a los TM 4,5,6..
               open(unit=21,file="matrices/RiaTM_Icampo"//trim(icampoch),form='unformatted', status = 'replace')
               write(21) RiaTM
               close(21)
               open(unit=21,file="matrices/RiaTMconj_Icampo"//trim(icampoch),form='unformatted', status = 'replace')
               write(21) RiaTM_conj
               close(21)
            End Do

         End do
      End do

      deallocate(Ri,R1,R2,vec1, vec2, vec3,vec4)
      deallocate(RiaTE, RiaTE_conj, RiaTM, RiaTM_conj )
      deallocate(block_sizeTE,block_sizeTM, offset_l_TE, offset_l_TM)


   end subroutine integral_radial
!!!!!!!!CORRECION CON EL ICAMPO DEBIDO A QUE EL LOOP SE HACE COMO L N M --> LA TUPLA SE ORDENA LEXICOGRAFICAMENTE COMO (L N M)
      !PROBLEMA PORQUE DESPUES EN BUILD MATRIX SE BARRE EN m EN MODES PERO ESA BASE ESTA ORDENADA COMO (L M N). ARREGLADO


   subroutine BuildMatrix()
      use complex_module
      use ModuleParameters
      use pi_module
      use sph_bessel_roots
      use ModuleState
      use clight_module
      implicit none
      integer :: i,j , l, n, m, imi, ima, jmi, jma, iostat, k
      real*8 :: b_a, b_adaga, wln, nte_ln, wln_p, ntm_ln , knl_p
      integer, allocatable :: vector_fock(:), estado_mat(:)
      integer, allocatable :: fock_izq(:), mat_izq(:),fock_der(:), mat_der(:)
      real*8 :: func

      ! complex*16, allocatable :: Hint(:,:)


      character(2)        :: icampoch

      character(len=8)    :: fmt1,fmt2
      character(len=32)   :: Rmaxch

      real*8, allocatable :: RiaTE(:,:), RiaTE_conj(:,:),RiaTM(:,:), RiaTM_conj(:,:)

      character(len=*), parameter :: fmt = '(2I4, 2(" (",F6.2,",",F6.2,")"))'

      
      integer :: lwork, lwork_opt, info, lda
      real*8, allocatable :: RWORK(:)
      complex*16, allocatable :: WORK(:)
      complex*16 :: WORK_QUERY(1)


      fmt1 = '(I3.3)'
      fmt2 = '(I2.2)'
 
   
      if(allocated(Hint)) deallocate(Hint)
      if(allocated(Hm)) deallocate(Hm)
      if(allocated(Hf)) deallocate(Hf)
      if(allocated(diag)) deallocate(diag)

      allocate(Hint(Nsystem,Nsystem), Hf(Nsystem), Hm(Nsystem), diag(Nsystem))

      Hint = c0
      Hm = 0.d0
      Hf = 0.d0
      diag = 0.d0

      write(*,*) "BuildTotalMatrix"

      ! call InitConfSpace()

      !H de fotones
      !Energia de una configuracion de fotones
      E0 = 0.d0

      Do i = 1, M_TE
         l = modes_matrix(i,2) !es el l en la matriz de estados fotonicos
         n = modes_matrix(i,4) !es el n en la matriz de estados fotonicos
         E0 = E0 + zBessTE(l, n)
      End do


      !falta definir dj_j_root para que encuentre las raices para los TM
      Do i = M_TE + 1, M_TE + M_TM
         l = modes_matrix(i,2) !es el l en la matriz de estados fotonicos
         n = modes_matrix(i,4) !es el n en la matriz de estados fotonicos
         E0 = E0 + zBessTM(l, n)
      End do

      E0 = E0*(1.d0/Rmax)*clight*0.5d0

      write(*,*)'E0 =', E0

      allocate(vector_fock(M_TE+M_TM),estado_mat(3))

      !!PARTE DIAGONAL
      Do i = 1, Nsystem

         call get_state(i,NsysMA,vector_fock,estado_mat) !aca ya tengo los estados fotonico y materiales

         Hf(i) = E_N(vector_fock) + E0 !construyo la matriz diagonal de las energias de los fotones (es un vector directamente)

         !lmn materiales
         l = estado_mat(1)
         n = estado_mat(3)
         Hm(i) = Energias(n,l) !!cada estado i-esimo viene con su Emat dada por get_state

         Hint(i,i) = c1*Hm(i) + c1*Hf(i) !!leno la parte diagonal en Hint
         diag(i)= Hm(i) + Hf(i)
      End do



      allocate(fock_izq(M_TE + M_TM), mat_izq(3))
      allocate(fock_der(M_TE + M_TM), mat_der(3))

      allocate(RiaTE(NsysMA,NsysMA),RiaTE_conj(NsysMA,NsysMA),RiaTM(NsysMA,NsysMA),RiaTM_conj(NsysMA,NsysMA))

      !H DE INTERACCION
      Do i = 1, NsysPH
         imi = (i-1)*NsysMA+1 !este es el inicio de las filas del bloque para esa conf de fock
         ima = i*NsysMA       !este es el final de las filas para ese bloque
         call get_state(imi, NsysMA,fock_izq, mat_izq )
         ! write(*,*)'fila',imi
         ! write(*,*)'estado fotonico de fock de izquerda(fila)', fock_izq  

         Do j = 1, NsysPH
            jmi = (j-1)*NsysMA+1 !este es el inicio de las columnas del bloque para esa conf de fock
            jma = j*NsysMA       !este es el inicio de las columanas para ese bloque
            call get_state(jmi, NsysMA,fock_der, mat_der)
            ! write(*,*)'columna',jmi
            ! write(*,*)'estado fotonico de fock de derecha(columna)', fock_der 

            Do m = 1, M_TE  !!este barre sobre modos fotonicos dentro de un dado bloque, en vez de checkear todos los modos, veo hasta el primero que no es nulo
               !checkeo si son cero o no para un dado m
               !obtengo el l correspondiente

               l = modes_matrix(m,2)
               n = modes_matrix(m,4)
               b_a = delta_a(fock_der, fock_izq, m) !destruccion
               b_adaga = delta_adaga(fock_der, fock_izq, m) !creacion


               write(icampoch,fmt2) m
               ! write(*,*) 'ESTADO m: ',modes_matrix(m,:)

               !frecuencias
               wln = zBessTE(l,n)*(1.d0/Rmax)*clight
               nte_ln = NTE(l,n) !!factor de normalizacion

               if(dabs(b_a).gt.0.000000001d0)then
                  open(unit=15, file="matrices/RiaTE_Icampo"//trim(icampoch), form='unformatted', action='read', iostat=iostat)
                  if (iostat /= 0) then
                     print *, "Error abriendo el archivo"
                     stop
                  endif
                  read(15) RiaTE
                  close(15)

                  Hint(imi:ima,jmi:jma) = Hint(imi:ima,jmi:jma) + lambda*sqrt(2.d0*pi/(wln))*(nte_ln)*(1.d0/sqrt(l*(l+1)*1.d0))*(-ci*(1.d0/(clight*me))*RiaTE*b_a ) !!aca pega esos bloques en la matriz grande en los rangos indicados
               ! exit !!esto es porque al primero distinto de cero, ya se que el resto sera cero

               endif


               if(dabs(b_adaga).gt.0.000000001d0)then
                  open(unit=15, file="matrices/RiaTEconj_Icampo"//trim(icampoch), form='unformatted', action='read', iostat=iostat)
                  if (iostat /= 0) then
                     print *, "Error abriendo el archivo"
                     stop
                  endif
                  read(15) RiaTE_conj
                  close(15)
                  Hint(imi:ima,jmi:jma) = Hint(imi:ima,jmi:jma) + lambda*sqrt(2.d0*pi/(wln))*(nte_ln)*(1.d0/sqrt(l*(l+1)*1.d0))*(- ci*(1/clight*me)*RiaTE_conj*b_adaga) !!aca pega esos bloques en la matriz grande en los rangos indicados
                  ! exit !!esto es porque al primero distinto de cero, ya se que el resto sera cero
               endif
               ! Hint(imi:ima,jmi:jma) = Hint(imi:ima,jmi:jma) + lambda*sqrt(1.d0/(2*wln))*(nte_ln)*(1.d0/sqrt(l*(l+1)*1.d0))*( -ci*(1.d0/(clight*me))*RiaTE*b_a - ci*(1/clight*me)*RiaTE_conj*b_adaga) !!aca pega esos bloques en la matriz grande en los rangos indicados
               ! exit !!esto es porque al primero distinto de cero, ya se que el resto sera cero


            End do

            !!!!!!!!!!!!!!!!!!! sigo la parte magnetica 
            Do m = M_TE+1,M_TE + M_TM  !!este barre sobre modos fotonicos dentro de un dado bloque, en vez de checkear todos los modos, veo hasta el primero que no es nulo
            !!!puedo hacer esto porque ya se que estan ordenado todos lo TE primeros y los TM despues. Es mas eficiente que poner if alpha = 1 o 2
               l = modes_matrix(m,2)
               n = modes_matrix(m,4)
               b_a = delta_a(fock_der, fock_izq, m) !destruccion
               b_adaga = delta_adaga(fock_der, fock_izq, m) !creacion


               write(icampoch,fmt2) m

   
               !frecuencias
               wln_p = zBessTM(l,n)*(1.d0/Rmax)*clight
               knl_p = zBessTM(l,n)*(1.d0/Rmax)
               ntm_ln = NTM(l,n) !!factor de normalizacion

               ! Hint(imi:ima,jmi:jma) = Hint(imi:ima,jmi:jma) + lambda*(1/knl_p)*sqrt(1.d0/(2*wln_p))*(ntm_ln)*(1.d0/sqrt(l*(l+1)*1.d0))*( ci*(1.d0/(clight*me))*RiaTM*b_a + ci*(1/clight*me)*RiaTM_conj*b_adaga) !!aca pega esos bloques en la matriz grande en los rangos indicados
               ! exit !!esto es porque al primero distinto de cero, ya se que el resto sera cero

               if(dabs(b_a).gt.0.000000001d0)then
                  open(unit=21, file="matrices/RiaTM_Icampo"//trim(icampoch), form='unformatted', action='read', iostat=iostat)
                  if (iostat /= 0) then
                     print *, "Error abriendo el archivo"
                     stop
                  endif
                  read(21) RiaTM
                  close(21)

                  Hint(imi:ima,jmi:jma) = Hint(imi:ima,jmi:jma) + lambda*sqrt(2.d0*pi/(wln))*(ntm_ln)*(1.d0/sqrt(l*(l+1)*1.d0))*(ci*(1.d0/(clight*me))*RiaTM*b_a ) !!aca pega esos bloques en la matriz grande en los rangos indicados
               ! exit !!esto es porque al primero distinto de cero, ya se que el resto sera cero

               endif


               if(dabs(b_adaga).gt.0.000000001d0)then
                  open(unit=21, file="matrices/RiaTMconj_Icampo"//trim(icampoch), form='unformatted', action='read', iostat=iostat)
                  if (iostat /= 0) then
                     print *, "Error abriendo el archivo"
                     stop
                  endif
                  read(21) RiaTM_conj
                  close(21)
                  Hint(imi:ima,jmi:jma) = Hint(imi:ima,jmi:jma) + lambda*sqrt(2.d0*pi/(wln))*(ntm_ln)*(1.d0/sqrt(l*(l+1)*1.d0))*( ci*(1/clight*me)*RiaTM_conj*b_adaga) !!aca pega esos bloques en la matriz grande en los rangos indicados
                  ! exit !!esto es porque al primero distinto de cero, ya se que el resto sera cero
               endif
               ! Hint(imi:ima,jmi:jma) = Hint(imi:ima,jmi:jma) + lambda*sqrt(1.d0/(2*wln))*(nte_ln)*(1.d0/sqrt(l*(l+1)*1.d0))*( -ci*(1.d0/(clight*me))*RiaTE*b_a - ci*(1/clight*me)*RiaTE_conj*b_adaga) !!aca pega esos bloques en la matriz grande en los rangos indicados
               ! exit !!esto es porque al primero distinto de cero, ya se que el resto sera cero



            End do

         End do
      End do
      
      deallocate(RiaTM,RiaTM_conj)

      !!! para escribir la matriz de Hint
      write(Rmaxch,fmt = '(F12.8)') Rmax
      Rmaxch = adjustl(Rmaxch)
      open(unit=17, file="resultados/matriz"//trim(Rmaxch)//".dat", status="replace", action="write")
      Do i = 1, Nsystem
         write(17,'( *(1x,ES25.16) )') (abs(Hint(i,j)), j=1,Nsystem)
      End Do
      close(17)


      !!!!! para checkear hermiticidad
      func = 0.d0
      k = 0
      Do i = 1, Nsystem
         Do j= i + 1, Nsystem
               if(abs(Hint(i,j)).gt.1.d-23)then
               func = func + abs(Hint(i,j)) + abs(Hint(j,i))
               k = k+1
               end if
         End Do
      End Do

      write(*,*) "orden de magnitud promedio no diagonal",func/(2.d0*k)

      func = 1.d-3*func
      Do i = 1, Nsystem
         Do j= i, Nsystem
               if(abs(Hint(i,j)-conjg(Hint(j,i))).gt.func)then
               write(*,*) "*******************************************"
               write(*, *) i,j
               write(*, *) Hint(i,j)
               write(*, *) Hint(j,i)
               write(*,*) "H no es hermítico. STOP"
               stop
               endif
         End Do
      End Do

   


      !DIAGONALIZACION DE MATRIZ HERMITICA!!!!!!!!!!!!!!!

      allocate(W(Nsystem))

      lda = Nsystem

      allocate(RWORK(max(1,3*Nsystem-2)))

      lwork = -1

      call zheev('V', 'U', Nsystem, Hint, lda, W, WORK_QUERY, lwork, RWORK, info)

      !lwork tiene el tamano optimo
      lwork_opt = int(real(WORK_QUERY(1)))

      allocate(WORK(lwork_opt))

      ! ---- Diagonalización real ----
      call zheev('V', 'U', Nsystem, Hint, lda, W, WORK, lwork_opt, RWORK, info)  !!--> aca estan ya los W autovalores en moduleParameters



      deallocate(WORK, RWORK)


   End subroutine BuildMatrix


   ! subroutine energia_conf_fock() !!esta es para guardar las energias de las configuraciones de fock que hay
   !    use moduleParameters
   !    use ModuleState

   !    implicit none

   !    integer :: i
   !    integer :: conf_fock(M_TE+M_TM)
   !    real*8 :: energia_fotones
      

   !    conf_fock = 0.d0

   !    if (allocated(EE)) deallocate(EE) 
   !    allocate(EE(NsysPH_TE*NsysPH_TM))
   !    EE=0.d0

   !    do i = 1, NsysPH_TE*NsysPH_TM
   !       conf_fock = state_matrix_ph(i,:)
   !       energia_fotones = E_N(conf_fock) + E0 !!depende de R
   !       EE(i) = energia_fotones 
   !    end do

   ! end subroutine energia_conf_fock




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! subroutine InitConfigurationSpace
   ! use ModuleParameters
   ! implicit none
   ! integer             :: ik, il, im, iMte, inte, iMtm, intm, isys
   !
   !
   ! allocate(kmat(1:Nsystem), lmat(1:Nsystem), mmat(1:Nsystem))
   ! allocate(Mte(1:Nsystem), nte(1:Nsystem), Mtm(1:Nsystem), ntm(1:Nsystem))
   !
   ! isys = 1
   !
   ! Do ik =  1, Nmax
   !    kmat(isys) = ik
   !    Do il =  0, Lmax
   !       lmat(isys) = il
   !       Do im =  -il, il
   !          mmat(isys) = im
   !          Do iMte = 1, Nmodes_TE
   !             Mte(isys) = iMte
   !             Do inte = 1, Nph_TE
   !                nte(isys) = inte
   !                Do iMtm = 1, Nmodes_TM
   !                   Mtm(isys) = iMtm
   !                   Do intm = 1, Nph_TM
   !                      ntm(isys) = intm
   !                      isys = isys + 1
   !                   End Do
   !                End Do
   !             End Do
   !          End Do
   !       End do
   !    End do
   ! End do
   ! End subroutine InitConfigurationSpace






















   ! subroutine InitCoulombWall
   ! use ModuleParameters
   ! implicit none
   ! integer             :: L, i, n
   ! real*8              ::  En,f0,f,fp,g,gp
   ! REAL*8              :: x1,x2
   ! real*8, allocatable :: x(:),w(:)
   !
   ! allocate(Energias(1:Nmax,0:Lmax))
   ! Energias = 0.d0
   ! Do L = 0, Lmax
   !    Do n = 1, Nmax
   !      Energias(n,L) =  -me*(Ztilde**2.d0)/(2.d0*n*n)
   !    enddo
   ! enddo
   !
   ! call NormContinuum

   ! endsubroutine



   ! subroutine NormContinuum
   ! use ModuleParameters
   ! implicit none
   ! integer             :: L, i, n
   ! real*8              ::  En,f0,f,fp,g,gp
   ! REAL*8              :: x1,x2
   ! real*8, allocatable :: x(:),w(:)
   !
   ! allocate(NormCoul(1:Nmax,0:Lmax))
   ! NormCoul    = 0.d0
   ! allocate(x(1:Nquad),w(1:Nquad))
   ! x1 = 0.d0
   ! x2 = Rmax
   ! call gauleg(x1,x2,x,w,Nquad)
   ! Do L = 0, Lmax
   !    Do n = 1, Nmax
   !       f0 = 0.d0
   !       En = Energias(n,L)
   !       Do i = 1, Nquad
   !          call scoul(Ztilde*sqrt(me),En,L,x(i)*sqrt(me),f,fp,g,gp,err)
   !          f0 = f0 + f*f*w(i)
   !       enddo
   !       NormCoul(n,L) = f0**(-0.5d0)
   !    enddo
   ! enddo
   ! deallocate(x,w)
   ! endsubroutine




   ! subroutine InitializeBessels
   ! use ModuleParameters
   ! implicit none
   ! integer             :: L, i, NM
   ! real*8              :: r, s1, sn
   ! real*8, allocatable :: Jbess(:), DJBess(:)
   ! real*8, allocatable :: Jlr(:,:), DJlr(:,:)
   ! real*8, allocatable :: a(:), b(:), c(:), d(:)
   !
   ! LmaxJ = Lmax + L0
   !
   ! allocate(Jbess(0:LmaxJ), DJBess(0:LmaxJ))
   ! allocate(Jlr(0:LmaxJ,1:Nrbound))
   !
   ! Jlr = 0.d0
   ! DJlr = 0.d0
   ! Jlr(0,1) = 1.d0
   ! do i = 2, Nrbound
   !    r = rsbound(i)
   !    call SPHJ(LmaxJ,k0*r,NM,Jbess,DJBess)
   !    Jlr(:,i)  = Jbess
   ! enddo
   ! deallocate(Jbess, DJBess)
   !
   ! allocate(a(1:Nrbound),b(1:Nrbound),c(1:Nrbound),d(1:Nrbound))
   ! allocate(jsa(0:LmaxJ,1:Nrbound),jsb(0:LmaxJ,1:Nrbound))
   ! allocate(jsc(0:LmaxJ,1:Nrbound),jsd(0:LmaxJ,1:Nrbound))
   ! ! allocate(djsa(0:LmaxJ,1:Nrbound),djsb(0:LmaxJ,1:Nrbound))
   ! ! allocate(djsc(0:LmaxJ,1:Nrbound),djsd(0:LmaxJ,1:Nrbound))
   !
   ! ! excepto si L=0 donde es 1 para la spherical
   ! ! SON LAS SECOND, AVERIGUAR, TIENEN QUE SER NULAS
   !
   !
   ! allocate(Jbess(1:Nrbound))
   !
   ! do L = 0, LmaxJ
   !       Jbess(:) = Jlr(L,:)
   !       s1 = 0.d0
   !       sn = 0.d0
   !       call spline(rsbound,Jbess,a,b,c,d,s1,sn,Nrbound)
   !       jsa(L,:) = a(:)
   !       jsb(L,:) = b(:)
   !       jsc(L,:) = c(:)
   !       jsd(L,:) = d(:)
   ! enddo
   !
   ! deallocate(Jbess)
   ! deallocate(Jlr)
   ! deallocate(a,b,c,d)
   ! End subroutine



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! subroutine InitializeBoundCoulomb
   ! use ModuleConstants
   ! use ModuleParameters
   ! implicit none
   ! integer, parameter :: ndim = 10000
   ! real*8 rad(ndim),p(ndim),q(ndim)
   ! integer ngp,ilast,ier
   ! common/radwf/rad,p,q,ngp,ilast,ier
   ! integer m, n
   ! real*8 r,f,g
   ! integer i,L
   ! real*8 En
   ! real*8 eps,dell
   ! real*8, allocatable :: a(:),b(:),c(:),d(:)
   ! real*8 s1,sn
   ! character(len=8) :: fmt ! format descriptor
   ! character(5) filenumber
   ! real*8,  external :: BoundCoulomb
   !
   ! Nrbound = ndim
   ! allocate(rsbound(1:Nrbound),rvsbound(1:Nrbound),Ebounds(1:Nbounds,0:Lmax))
   ! do i=1,Nrbound
   !   rsbound(i) = (Rmax*(i-1))/(1.d0*(Nrbound-1))
   ! enddo
   ! rvsbound = Z
   ! rad      = rsbound
   ! ngp      = Nrbound
   ! call vint(rsbound,rvsbound,Nrbound)
   ! Ebounds = 0.d0
   ! fmt = '(I3.3)'
   ! eps   = 0.001d0
   ! dell  = eps
   ! s1 = 0.d0
   ! sn = 0.d0
   ! allocate(a(1:Nrbound),b(1:Nrbound),c(1:Nrbound),d(1:Nrbound))
   ! allocate(bsa(1:Nbounds,0:Lmax,1:Nrbound),bsb(1:Nbounds,0:Lmax,1:Nrbound))
   ! allocate(bsc(1:Nbounds,0:Lmax,1:Nrbound),bsd(1:Nbounds,0:Lmax,1:Nrbound))
   ! do L = 0, Lmax
   !    do n = L + 1, Nbounds
   !       !write(*,*) "L",L,"n",n
   !       En = -Z*Z/(2.d0*mu*((n)**2))
   !       call sbound(En,eps,dell,n,L)
   !       if(ier.ne.0)then
   !          write(*,*) "sbound, ier",ier
   !       endif
   !       Ebounds(n,L) = En
   !       call spline(rad,p,a,b,c,d,s1,sn,Nrbound)
   !       bsa(n,L,:) = a(:)
   !       bsb(n,L,:) = b(:)
   !       bsc(n,L,:) = c(:)
   !       bsd(n,L,:) = d(:)
   !    enddo
   ! enddo
   ! deallocate(rvsbound,a,b,c,d)
   ! End subroutine InitializeBoundCoulomb

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! subroutine findZerosCoulombBox
   ! use ModuleConstants
   ! use ModuleParameters
   ! implicit none
   ! integer L
   ! real*8  E0,E1,E2,DeltaE,fthresh,f0,f1,f2,fp,g,gp
   ! integer i
   !
   ! fthresh = 10.d0**(-7.d0)
   !
   ! allocate(Energias(1:Nmax,0:Lmax))
   ! Energias = 0.d0
   !
   ! Do L = 0, Lmax
   !    E0      = pi*pi/(2.d0*mu*Rmax)  !Energía fundamental en la caja de tamaño Rmax
   !    E0      = E0/10.d0
   !    DeltaE  = E0
   !    i = 1
   !    call scoul(Z*sqrt(mu),DeltaE,L,Rmax*sqrt(mu),f0,fp,g,gp,err)
   !
   !    if (abs(f0).le.fthresh) then
   !       Energias(i,L) = DeltaE
   !       i = i+1
   !    endif
   !
   !    E0 = E0 + DeltaE
   !    E1 = E0 + DeltaE
   !
   !    call scoul(Z*sqrt(mu),E1,L,Rmax*sqrt(mu),f1,fp,g,gp,err)
   !    1  continue
   !    call scoul(Z*sqrt(mu),E0,L,Rmax*sqrt(mu),f0,fp,g,gp,err)
   !    11  continue
   !
   !    if(f0*f1 < 0.d0)then
   !       2  continue
   !       E2 = E0+((E1-E0)/2.d0)
   !       call scoul(Z*sqrt(mu),E2,L,Rmax*sqrt(mu),f2,fp,g,gp,err)
   !       if (abs(f2).le.fthresh) then
   !          !write(*,*) "encontro",L,i,E0
   !          !pause
   !          Energias(i,L) = E0
   !          if (i.eq.Nmax) goto 3
   !          i = i + 1
   !          E0 = E2 + DeltaE
   !          call scoul(Z*sqrt(mu),E0,L,Rmax*sqrt(mu),f0,fp,g,gp,err)
   !          E1 = E0 + DeltaE
   !          call scoul(Z*sqrt(mu),E1,L,Rmax*sqrt(mu),f1,fp,g,gp,err)
   !          goto 11
   !       endif
   !
   !       if(f2*f1 < 0.d0)then
   !          E0 = E2
   !          f0 = f2
   !          goto 2
   !       else
   !          f1 = f2
   !          E1 = E2
   !          goto 2
   !       endif
   !    else
   !
   !       E0 = E1
   !       f0 = f1
   !       E1 = E0 + DeltaE
   !       call scoul(Z*sqrt(mu),E1,L,Rmax*sqrt(mu),f1,fp,g,gp,err)
   !       goto 11
   !    endif
   !
   !    3  continue
   ! End Do
   ! endsubroutine



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! subroutine SaveBasisData
   ! use ModuleConstants
   ! use ModuleParameters
   ! implicit none
   ! logical :: dir_e
   ! inquire(file='./matrices/test', exist=dir_e)
   ! if (dir_e) then
   ! !   write(*,*) "dir exists!",dir_e
   ! else
   !   ! workaround: it calls an extern program...
   ! !   write(*,*) "dir does not exists!",dir_e
   !   call system('mkdir -p ./matrices')
   !   Open(unit=10,file="matrices/test")
   !   close(10)
   ! end if
   ! open(unit=10,file="matrices/AlltheData",form='unformatted')
   ! write(10) Rmax
   ! write(10) R0max
   ! write(10) Z
   ! write(10) mu
   ! ! write(10) k0
   ! write(10) theta0
   ! write(10) varphi0
   ! write(10) N0
   ! write(10) L0
   ! write(10) M0
   ! write(10) Nmax
   ! write(10) Lmax
   ! write(10) Nrbound
   ! write(10) Nbounds
   ! write(10) err
   ! write(10) Nquad
   ! write(10) NormCoul
   ! write(10) Energias
   ! write(10) Ebounds
   ! write(10) rsbound
   ! write(10) bsa
   ! write(10) bsb
   ! write(10) bsc
   ! write(10) bsd
   ! ! write(10) jsa
   ! ! write(10) jsb
   ! ! write(10) jsc
   ! ! write(10) jsd
   ! close(10)
   ! endsubroutine
   !
   !
   ! subroutine ReadBasisData
   ! use ModuleConstants
   ! use ModuleParameters
   ! implicit none
   ! logical :: dir_e
   ! inquire(file='./matrices/test', exist=dir_e)
   ! if (dir_e) then
   ! !   write(*,*) "dir exists!",dir_e
   ! else
   !   ! workaround: it calls an extern program...
   ! !   write(*,*) "dir does not exists!",dir_e
   !   call system('mkdir -p ./matrices')
   !   Open(unit=10,file="matrices/test")
   !   close(10)
   ! end if
   ! open(unit=10,file="matrices/AlltheData",form='unformatted')
   ! read(10) Rmax
   ! read(10) R0max
   ! read(10) Z
   ! read(10) mu
   ! ! read(10) k0
   ! read(10) theta0
   ! read(10) varphi0
   ! read(10) N0
   ! read(10) L0
   ! read(10) M0
   ! read(10) Nmax
   ! read(10) Lmax
   ! read(10) Nrbound
   ! read(10) Nbounds
   ! allocate(Energias(1:Nmax,0:Lmax),NormCoul(1:Nmax,0:Lmax))
   ! allocate(rsbound(1:Nrbound),Ebounds(1:Nbounds,0:Lmax))
   ! allocate(bsa(1:Nbounds,0:Lmax,1:Nrbound),bsb(1:Nbounds,0:Lmax,1:Nrbound))
   ! allocate(bsc(1:Nbounds,0:Lmax,1:Nrbound),bsd(1:Nbounds,0:Lmax,1:Nrbound))
   ! ! LmaxJ = Lmax + L0
   ! ! allocate(jsa(0:LmaxJ,1:Nrbound),jsb(0:LmaxJ,1:Nrbound))
   ! ! allocate(jsc(0:LmaxJ,1:Nrbound),jsd(0:LmaxJ,1:Nrbound))
   !
   ! read(10) err
   ! read(10) Nquad
   ! read(10) NormCoul
   ! read(10) Energias
   ! read(10) Ebounds
   ! read(10) rsbound
   ! read(10) bsa
   ! read(10) bsb
   ! read(10) bsc
   ! read(10) bsd
   ! ! read(10) jsa
   ! ! read(10) jsb
   ! ! read(10) jsc
   ! ! read(10) jsd
   ! close(10)
   ! endsubroutine


end module subrutinas
