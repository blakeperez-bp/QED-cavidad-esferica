module spherical_bessel
  implicit none
contains
  function sph_bessel_j(l, x) result(val)
    implicit none
    integer, intent(in) :: l
    real(8), intent(in) :: x
    real(8) :: val
    real(8) :: j0, j1, jm1, jcurr
    integer :: n
  !!relacion de recurrencia del NIST  (https://dlmf.nist.gov/10.51)
    if (l == 0) then
       val = sin(x)/x
    else if (l == 1) then
       val = sin(x)/x**2 - cos(x)/x
    else
       j0 = sin(x)/x
       j1 = sin(x)/x**2 - cos(x)/x
       do n = 1, l-1
          jm1 = j0
          jcurr = j1
          j1 = (2*n+1)/x * jcurr - jm1
          j0 = jcurr
       end do
       val = j1
    end if
  end function sph_bessel_j


   real(8) function dj_rjl(l, x)   !!deriavada de r*j_l
      implicit none
      integer, intent(in) :: l
      real(8), intent(in) :: x
      real(8) :: jl, jlm1

      if (x == 0.d0) then
         ! comportamiento cerca de 0: d/dx[x j_l] ~ 0 para l>0, para l=0 es cos(0)=1
         if (l == 0) then
            dj_rjl = 1.d0
         else
            dj_rjl = 0.d0
         end if
         return
      end if

      jl = sph_bessel_j(l, x)
      if (l == 0) then
         ! usamos j_{-1} solo si lo definimos; mejor usar la forma alternativa:
         ! para l=0 conocemos la derivada exacta: d(x j0)/dx = cos x
         dj_rjl = cos(x)
      else
         jlm1 = sph_bessel_j(l-1, x)
         dj_rjl = x * jlm1 - l*1.d0* jl !!relacion dada por el Nist
      end if
   end function dj_rjl


   real*8 function dj_l(l,x) result(val)  !!derivada de la bessel relacion sacada del NIST
      implicit none
      integer, intent(in) :: l
      real*8, intent(in)  :: x

      if (l==0) then
         val = (x*cos(x) - sin(x))/x**2
      elseif ( l>0) then
         val = sph_bessel_j(l-1, x) - ((l+1)/x)*sph_bessel_j(l, x)
      end if

   end function


   real*8 function ddxj(l,x) result(val) !!derivada segunda de xj_l
      implicit none
      integer, intent(in) :: l
      real*8, intent(in)  :: x

      if (l==0) then
        val = -sin(x)
      elseif(l>0) then
         val = 2*dj_l(l,x) + x*( dj_l(l-1,x)  + ((l+1)/x**2)*sph_bessel_j(l,x) - ((l+1)/x)*dj_l(l,x)   )
      end if
   end function 


end module spherical_bessel
!!!!!todas estan funciones ya estan checkeadas que estan bien


!!revisar si esta bien el metodo, verlo en una tabla
module sph_bessel_roots
   use spherical_bessel
   use pi_module
   implicit none
   private
   public :: bessel_j_root
   public :: bessel_jprime_root

   contains

         function bessel_j_root(l, n) result(root)
            ! Devuelve la n-ésima raíz de j_l(x)
            integer, intent(in) :: l, n
            real(8) :: root
            real(8) :: a, b, fa, fb, c, fc
            real(8), parameter :: tol = 1.0d-12
            integer :: max_iter, iter

            max_iter = 1000

            ! Estimación inicial usando aproximación asintótica
            a = (n-1)*acos(-1.d0) + l*acos(-1.d0)/2.d0
            b = n*acos(-1.d0) + l*acos(-1.d0)/2.d0

            fa = sph_bessel_j(l,a)
            fb = sph_bessel_j(l,b)

            if (fa*fb > 0.d0) then
               print *, "Error: no hay cambio de signo en el intervalo!"
               root = -1.d0
               return
            end if

            ! Bisección
            iter = 0
            do while (abs(b-a) > tol .and. iter < max_iter)
               c = (a+b)/2.d0
               fc = sph_bessel_j(l,c)
               if (fa*fc < 0.d0) then
                  b = c
                  fb = fc
               else
                  a = c
                  fa = fc
               end if
               iter = iter + 1
            end do

            root = (a+b)/2.d0

         end function bessel_j_root

       
         function bessel_jprime_root(l, n) result(root) !!para sacara las raices para las frecuencias de TM
            implicit none
            integer, intent(in) :: l, n
            real(8) :: root
            real(8), parameter :: tol = 1.0d-12
            integer, parameter :: max_iter = 2000
            integer :: max_samples, i, iter
            real(8) :: x1, x2, dx, a, b, fa, fb, c, fc
            logical :: found

            if (n < 1) then
               print *, "bessel_jprime_root: n debe ser >= 1"
               root = -1.d0
               return
            end if

            if (l == 0) then
               ! caso exacto: cos x = 0 -> x = (2n-1)*pi/2
               root = dble(2*n - 1) * pi / 2.d0
               return
            end if

            ! Inicio de la búsqueda: estimación asintótica para ceros puede usarse como guía,
            ! pero aquí muestreamos desde un valor pequeño hacia la derecha.
            x1 = 1.d-6
            dx = pi/4.d0
            x2 = x1 + dx
            fa = dj_rjl(l, x1)
            fb = dj_rjl(l, x2)

            found = .false.
            max_samples = 20000
            i = 0
            do while (.not. found .and. i < max_samples)
               if (fa*fb <= 0.d0) then
                  a = x1
                  b = x2
                  found = .true.
                  exit
               end if
               x1 = x2
               x2 = x2 + dx
               fa = fb
               fb = dj_rjl(l, x2)
               i = i + 1
            end do

            if (.not. found) then
               print *, "bessel_jprime_root: no se encontro cambio de signo para l=", l, " n=", n
               root = -1.d0
               return
            end if

            ! Si queremos la n-ésima raíz, repetimos el proceso hasta que encontremos n cambios de signo.
            ! (Aquí el muestreo horizontal encontró el primer intervalo—para la n-ésima raíz, 
            !  habría que seguir contando intervalos con cambio de signo.)
            ! Implementación: contar cambios de signo hasta n
            ! (versión simple: reiniciamos la búsqueda y contamos)
            x1 = 1.d-6
            dx = pi/4.d0
            fa = dj_rjl(l, x1)
            i = 0
            do while (i < n)
               x2 = x1 + dx
               fb = dj_rjl(l, x2)
               if (fa*fb <= 0.d0) then
                  ! Bisección en [x1,x2] para localizar la raíz
                  a = x1; b = x2
                  iter = 0
                  do while (abs(b - a) > tol .and. iter < max_iter)
                     c = 0.5d0*(a + b)
                     fc = dj_rjl(l, c)
                     if (fa*fc <= 0.d0) then
                     b = c
                     fb = fc
                     else
                     a = c
                     fa = fc
                     end if
                     iter = iter + 1
                  end do
                  i = i + 1
                  if (i == n) then
                     root = 0.5d0*(a + b)
                     return
                  end if
                  ! continuar buscando hacia la derecha
                  x1 = x2
                  fa = fb
               else
                  x1 = x2
                  fa = fb
               end if
            end do

            ! Si salimos sin encontrar:
            root = -1.d0
         end function bessel_jprime_root


   end module sph_bessel_roots








