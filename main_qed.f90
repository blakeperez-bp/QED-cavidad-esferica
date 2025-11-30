program SphericalQED
! use complex_module
use ModuleParameters
use funciones_index
use ModuleState
use funciones_angulares
use sph_bessel_roots
use subrutinas
use clight_module
use checkk

!use MKL_DFTI
implicit none

!declaro variables
integer     i, j, n, m, l
real(8)  :: t, t1, t2
character(8)     :: timeCh
real*8           :: func0
real(8) :: a, b, z, result
integer :: ierr, k


integer, allocatable :: conf_fock(:), est_material(:)

real*8 :: II3_, II4



integer :: der(4),izq(4)
real*8 :: bb, b_d


integer, parameter:: p = 1 !!cantidad de pasos para el vector de R's
real*8 :: start_val, end_val
real*8 :: R(p)

real*8, allocatable :: all_eigs(:,:)
real*8, allocatable :: mat(:,:), fotones(:,:), diagonal(:,:)



!real*8, external :: MargPsiX, MargLambdaT,MargLambdaTLARGE0, MargLambdaTLARGE
write(*,*) "SphericalQED"


SavedONOFF = "newcalc"
!SavedONOFF = "openold"


if(SavedONOFF.eq."newcalc")then
!set up parameters

   !Radio de la cavidad, Número máximo de funciones radiales, momento angular máximo de funciones materiales

   Nmax       = 1    ! <  Nfbox:  Número de funciones radiales materiales por cada  0 < = l  < = Lmax
   Lmax       = 1     !(2l+1) funciones por l, dando un total de (1+Lmax)^{2}
   Nfbox      = 100   ! Numero de funciones de base analíticas sin(n*pi*r/R)*(2/R)^{1/2}  con las que vo a representar a mis estados radiales
   if(Nmax.gt.Nfbox)then
      write(*,*) "Nmax debe ser menor o igual a Nfbox"
      stop
   EndIf
   me         = 1.d0   !en caso que haga falta vestirla (no implementado)
   Ze         = -1.d0
   Zn         = 3.d0
   Ztilde     = Ze/4.d0

   ! Numero de modos transversal eléctrico y magnético, y numero máximo de fotones en cada uno

   Rmodes_TE  = 1  ! numero de componentes radiales
   Lmodes_TE  = 1  ! número de valores posibles de 1 \le ele \le LmodesMax: componentes en ele, cada una de ellas con 2l+1 emes
   Nph_TE     = 1  ! número máximo de fotones en cada modo

   Rmodes_TM  = 1  ! numero de componentes radiales
   Lmodes_TM  = 1  ! número de valores posibles de 1 \le ele \le LmodesMax: componentes en ele, cada una de ellas con 2l+1 emes
   Nph_TM     = 1  ! número máximo de fotones en cada modo

   m0 = 0         ! -min{l,m0} <= m <= min{l,m0} para modos fotonicos

   lambda = 1.d0 !!-->factor de escaleo para el acoplamiento
!    lambda = 0.d0 !!-->factor de escaleo para el acoplamiento

   !Tamaño del sistema
   ! M_TE = Rmodes_TE*Lmodes_TE*(2 + Lmodes_TE) !!cantidad de modos  electricos
   ! M_TM = Rmodes_TM*Lmodes_TM*(2 + Lmodes_TM) !!cantidad de modos magneticos
   
   ! CON m0
   M_TE = ((2*m0+1)*Lmodes_TE - m0*(m0-1))*Rmodes_TE
   M_TM = ((2*m0+1)*Lmodes_TM - m0*(m0-1))*Rmodes_TM

   NsysPH_TE = (Nph_TE+1)**M_TE
   NsysPH_TM = (Nph_TM+1)**M_TM

   write(*,*)'NsysPH_TE, NsysPH_TM ', NsysPH_TE, ' ', NsysPH_TM

   write(*,*) "Cantidad de modos de campo:"
   write(*,*) "de TE: ", M_TE, " de TM: ", M_TM
   NsysPH    =          (Nph_TE + 1)**(M_TE)
   NsysPH    = NsysPH*((Nph_TM + 1)**(M_TM))
   write(*,*) "Configuraciones de Fock:", NsysPH
   NsysMA    = Nmax*(1+Lmax)*(1+Lmax)
   write(*,*) "Modos materiales:", NsysMA
   Nsystem    = NsysMA*NsysPH
   write(*,*) "Total de configuraciones:",Nsystem
   Nquad      = 2000 !Numero de puntos de cuadratura


   
   call InitZerosBessel()
   call int_auxiliar_material() !hay que hacerlo una vez sola, a menos que se cambie Nfbox o Lmodes/Rmodes --> esto crea las matrices

   call leer_matrices()   !--> esto es para leerlas y dejarlas guardadas en ModuleParameters

   call InitConfSpace()  !!inicializo el espacio de estados materiales x fotonicos

   
   call angular_archivo() !!esta me crea y guarda un vector ordenado con todos los factores angulares, para una dada conf de parametros se hace solo una vez
   call leer_angular() !! este lo lee y lo guarda en module parameters para que ya quede disponible




   ! start_val = 50.0d0
   ! end_val = 100.0d0

   ! R = [(start_val + (end_val - start_val) * real(i-1) / real(p-1), i = 1, p)]
   R = [50.d0]



!    open(unit=16, file="resultados/autovalores.dat", status="replace", action="write") !!abrimos el archivo para general la matriz de autoenergias, cada columna es un Rmax
   ! open(unit=16, file="resultados/autovalores.dat", status="unknown", action="write", position="append")
   ! open(unit=18, file="resultados/energia_mat.dat", status="replace", action="write")  !!son las energias de la parte material
   ! open(unit=19, file="resultados/energia_fotones.dat", status="replace", action="write") !!son las energias de las conf de fock
   ! open(unit=84, file="resultados/diagonal.dat", status="unknown", action="write", position="append")
!    open(unit=84, file="resultados/diagonal.dat", status="replace", action="write")

   allocate(all_eigs(Nsystem, size(R) )) !!matriz que guarda en sus columnas a los autovalores para dado R
   ! allocate(mat(Nsystem, size(R)))  !!para las energias de la parte material solamente
   ! allocate(fotones(NsysPH_TE*NsysPH_TM, size(R)))

   allocate(diagonal(Nsystem, size(R) ))


   do i = 1, size(R)
      Rmax = R(i)

      write(*,*)'VALORE DE RMAX', Rmax,i, size(R)

      call ctes_normalizacion()


      call InitMaterialFunctions()

      call integral_radial()

      call BuildMatrix()


      all_eigs(:,i) = W(:)
      diagonal(:,i) = diag(:)

!       write(16,*) Rmax, (W(j), j=1,10)
      ! Útil si los valores tienen rangos muy diferentes
      ! write(16,'(F12.6, 64ES16.8)') Rmax, (W(j), j=1,64)
      ! write(84,'(F12.6, 64ES16.8)') Rmax, (diag(j), j=1,64)

      ! write(16,'(F20.10,1x, *(ES25.16,1x))') Rmax, (W(j), j=1,size(W))
      ! write(84,'(F20.10,1x, *(ES25.16,1x))') Rmax, (diag(j), j=1,size(diag))


      ! mat(:,i) = Hm(:)

      ! call energia_conf_fock()
      ! fotones(:,i) = EE(:)


      if (allocated(NTE) .and. allocated(NTM)) deallocate(NTE, NTM)
      if(allocated(W)) deallocate(W)
      ! if(allocated(Hm)) deallocate(Hm)
      if (allocated(diag)) deallocate(diag)
      if (allocated(Hint)) deallocate(Hint)

   end do


   close(16)
   close(84)

   ! do i = 1, Nsystem
   !    write(18,'( *(f15.6,1x) )') mat(i,:)
   ! end do

   ! close(18)


   ! do i = 1,  NsysPH_TE*NsysPH_TM
   !    write(19,'( *(f15.6,1x) )') fotones(i,:)
   ! end do

   ! close(19)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Open a file to store all energies
! open(100, file='energies_table_nucleo_cond.dat', status='replace')

! Write header
! write(100, '(A)', advance='no') 'Rmax    '
! do l = 0, Lmax
!   do n = 1, Nmax
!     write(100, '(A, I2, A, I1, A)', advance='no') '  E(', n, ',', l, ')'  ! E(10,0)
!   end do
! end do
! write(100, *)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write data
! do Rmax = 100.0d0,100.0d0

!    call InitMaterialFunctions
!    print *, "First energy value:", Energias(1,0)
!    write(100, '(F10.2)', advance='no') Rmax
!    do l = 0, Lmax
!       do n = 1, Nmax
!          write(100, '(F12.6)', advance='no') Energias(n, l)
!          write(*, '(F12.6)', advance='no') Energias(n, l)
!       end do
!    end do
!    write(100, *)
!    !deallocate(Energias, Estados, Hmaterial, En)  ! Clean up for next iteration
! end do
! deallocate(Energias)
!    pause


elseif(SavedONOFF.eq."openold")then
!    call Readdata
   !call InitGreenXEF
endif




End program SphericalQED


