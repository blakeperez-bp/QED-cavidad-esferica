

module complex_module
complex*16, parameter :: c1=(1d0,0d0)
complex*16, parameter :: ci=(0d0,1d0)
complex*16, parameter :: c0=(0d0,0d0)
real*8, parameter :: sq2i=0.70710678118654752440d0 !inversa de la raiz de dos
end module

module pi_module
real*8, parameter :: pi=3.1415926535897932d0
end module

module clight_module
    real*8, parameter :: alfa = 7.297352569e-3  !!extraido del PDG particle data group
    real*8, parameter :: clight = 1.d0/alfa
end module