      MODULE  CONSTANT
        USE PREC
        implicit none 
        real(dp), parameter :: b2a = 0.52917721092d0
        real(dp), parameter :: pi = 4.d0*DATAN(1.d0)
        real(dp), parameter :: invsqrt2pi = 1.d0/sqrt(2.0d0*pi)
        real(dp), parameter :: invsqrtpi = 1.d0/sqrt(pi)
        real(dp), parameter :: invpi = 1.d0/pi
        real(dp), parameter :: ev2h = 1.d0/27.21138505d0

      END MODULE
