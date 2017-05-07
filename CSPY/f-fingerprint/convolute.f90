 module convolute 
 use prec
 use constant
 contains 
  
  subroutine gaussian_smearing(x0, y0, xdata, ydata, npoints, sigma)
    !---------------------------------------------------------------------------
    !
    ! gaussian_smearing: a conventional convolution where the gaussian function
    ! of width sigma is the second factor in the integral.
    ! G(x) = 1/(sqrt(2*pi)*sigma) * exp(- x**2/(2*sigma**2)) 
    !---------------------------------------------------------------------------
    
    integer, intent(in) :: npoints
    real(dp), intent(in) :: x0, y0, sigma
    real(dp), intent(out),dimension(npoints) :: xdata, ydata
    real(dp), dimension(npoints) :: tmp
    integer :: i, j
    
    do i=1,npoints
        tmp(i) =  y0 * invsqrt2pi / sigma * dexp(      &
&                -0.5d0 * ((xdata(i) - x0) / sigma)**2)
    end do
    
    ydata = tmp
    
    return
  end subroutine
  
  subroutine lorentzian_convolute(x0, y0, xdata, ydata, npoints, sigma)
    !---------------------------------------------------------------------------
    !
    ! lorentzian_convolute: a conventional convolution where the lorentzian function
    ! of width sigma is the second factor in the integral.
    ! F(x) = 1/pi * (0.5*sigma/((x-x0)**2+ (0.5*sigma)**2)
    !---------------------------------------------------------------------------
    
    integer, intent(in) :: npoints
    real(dp), intent(in) :: x0, y0, sigma
    real(dp), intent(out),dimension(npoints) :: xdata, ydata
    real(dp), dimension(npoints):: tmp
    integer :: i, j
    
    do i=1,npoints
        tmp(i) =  y0 * invpi * (0.5d0*sigma / ((xdata(i) - x0)**2 + &
 &                (0.5d0*sigma)**2))
    end do
    
    ydata = tmp
    
    return
  end subroutine
  
  end module 
