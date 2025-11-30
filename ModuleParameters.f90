module ModuleParameters
!!aca estan todas las declaraciones de las variables


real*8, parameter :: err = 0.00000001d0
real*8  :: Ze, me, Ztilde, Zn
real*8  :: Rmax
integer :: Nsystem, NsysMA, NsysPH  !Tamaño del sistema
integer :: Nfbox, Nmax, Lmax   !(2l+1) funciones por l, dando un total de (1+Lmax)^{2}
! Numero de modos transversal eléctrico y magnético, y numero máximo de fotones en cada uno
integer ::    Rmodes_TE  ! numero de componentes radiales
integer ::    Lmodes_TE  ! número de valores posibles de 1 \le ele \le LmodesMax: componentes en ele, cada una de ellas con 2l+1 emes
integer ::    Nph_TE     ! número máximo de fotones en cada modo

integer ::     NsysPH_TE, NsysPH_TM

integer ::    Rmodes_TM  ! numero de componentes radiales
integer ::    Lmodes_TM  ! número de valores posibles de 1 \le ele \le LmodesMax: componentes en ele, cada una de ellas con 2l+1 emes
integer ::    Nph_TM     ! número máximo de fotones en cada modo
integer :: m0  !!valor que determina los valores de m para los modos fotonicos

integer :: M_TE, M_TM !!numero de modos electricos y magneticos

! matrices del sistema de autovalores material
real*8, allocatable :: Hmaterial(:,:)

!                              nrad, ele,  indice radial del seno(n*pi*r/R)
real*8, allocatable :: Estados(:,:,:)
real*8              :: Etilde
real*8 :: E0 !energia de vacio de fotones
real*8, allocatable :: En(:)
real*8, allocatable :: Energias(:,:)
real*8, allocatable :: zBessTE(:,:)
real*8, allocatable :: zBessTM(:,:)
real*8, allocatable :: NTE(:,:)
real*8, allocatable :: NTM(:,:)

real*8, allocatable :: EE(:)
real*8, allocatable :: Hm(:), Hf(:), diag(:) !--> Hamiltonianos

complex*16, allocatable :: Hint(:,:)

real*8,  allocatable :: W(:) !--> aca van a estar los autovalores de Hint para cada Rmax

real*8 :: lambda !--> factor de acoplamiento 


!quadratura
integer             :: Nquad
real*8, allocatable :: Txtom2(:,:), Tsurface(:,:), Tcoul(:,:), S(:,:), S1(:,:), S2(:,:)
real*8, allocatable :: Xquad(:), Wquad(:)

!angular
real*8 , allocatable:: factor_angular1(:), factor_angular2(:), factor_angular3(:), factor_angular4(:) 

! real*8, allocatable :: EE(:)


!norma, grids
character(len=7)    :: SavedONOFF

! Indices del sistema           k, l, m, Mte, nte , Mtm, ntm
integer, allocatable :: GlobalI(:,:,:,:,:), kmat(:), lmat(:), mmat(:), Mte(:), Mtm(:)

real*8, allocatable :: RiPlus(:,:), RiMinus(:,:)


! !integer :: ind, nvar
! integer, allocatable :: ranges (:,:)
! integer,allocatable:: range_ph (:,:)
! !integer, allocatable :: xs(:)
!
! !Todos estos de tamaño Nsystem
! integer, allocatable :: nMA(:), lMA(:), mMA(:)
! integer, allocatable :: npPH_TE(:), lPH_TE(:), mPH_TE(:), nRPH_TE(:)
! integer, allocatable :: npPH_TM(:), lPH_TM(:), mPH_TM(:), nRPH_TM(:)



! Coulombianas de la caja esférica
!real*8, allocatable :: NormCoul(:,:)

! real*8, allocatable :: GreenX(:,:), GreenT(:,:), GreenXT(:,:)
! real*8, allocatable :: MatrixXT(:,:), OverlapXT(:,:)
! real*8, allocatable :: EigenGreenX(:), EigenGreenT(:), EigenGreenXT(:)
! real*8, allocatable :: phi0XT(:), p0X(:), p0T(:)
! real*8, allocatable :: AutoX(:) , AutoT(:)
! real*8, allocatable :: NormXT(:), AutoXT(:)


!Index for the global time iteration
!integer        :: iTiter
!character(25)  :: Chtfilen

!Parameters for the basis and expansion of the jump probability hystogram fit
!integer :: NmaxH, MmaxH
!real*8 :: xminH, xmaxH, xmidH, xboxH, tminH, tmaxH, omegaOA
!real*8, allocatable :: Coefv(:,:)

!splines for x coordinate, for the sotution at tf=Tmax. The purpose is to generate a new initial condition for
!a calculation for t>Tmax
!real*8, allocatable :: tfbsa(:), tfbsb(:), tfbdc(:), tfbsd(:), rfi(:), ptfi(:)
!integer             :: Nbstf

!Solve using finite difference apprach. Whe store p for a runtime value t as a grid in x as PtithX(odis,Nquad),
!where the first argument odis stand for the sucessive t values we have to store for time discretization (for example odis=3 for a tree term recurrence
!obtained by a second order differentiation of the first derivative, odis>3 for runge kutta schemes) and Nquad stands for quadrature values in the x coordinate.
!
!integer             :: odis
!real*8              :: dt
!real*8, allocatable :: PtithX(:,:)

!real*8, allocatable :: IntdPsidtPxt(:) !saved integral from 0 to t, whose value is actualized by an amount dt and used in the f function for the RK4 method

!redefinition of the basis to a general-non orthonormal set. Normalization vectors
!real*8, allocatable :: NormX(:), NormT(:)
!real*8, allocatable :: OverX(:,:), OverT(:,:)

!splines for the time integral of the sacade time distribution
!real*8, allocatable :: tbsa(:), tbsb(:), tbsc(:), tbsd(:), tpsi(:), psi(:)


!real*8, allocatable :: jsa(:,:)  ,jsb(:,:)  ,jsc(:,:)  ,jsd(:,:)
!real*8, allocatable :: djsa(:,:) ,djsb(:,:) ,djsc(:,:) ,djsd(:,:)

!coefficients
!real*8,     allocatable :: Rcn_f(:,:,:), Rcn_b(:,:,:)
!complex*16, allocatable :: Cn_f(:,:,:), Cn_b(:,:,:)
!plotting
!integer                 :: NXplot, NTplot
!real*8                  :: Rplot , Zplot

   !Parametros de las distribuciones Marginales

!real*8                  :: varepsilon, x0, alpha, b0, x_0_deltaprox
!integer                 :: m0

end module



