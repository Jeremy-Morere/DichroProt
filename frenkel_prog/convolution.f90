subroutine convolution

!From the list of rotational stregth and transition energy
!convolute a gaussian or lorentzienne function

use declare
implicit none


npeak = size(energy_c)

if (type .eq. 'G') then
   c = fwhm
   q = 1.0d0/(2.294e1*pi)

   do i = 1,npeak
    if ((energy_c(i) .gt. borneinf) .and. (energy_c(i) .lt. bornesup)) then
      do j = 1,npoint
          absorb(j) = absorb(j) + rot_c(i)*dexp(-((lambda_c(j)-energy_c(i))/c)**2)
      enddo
    endif
   enddo
   absorb = absorb*q 

else if (type .eq. 'L') then
   c = fwhm
   q = 1.0d0/(2.294e1*pi) !cf Scott 2021

   do i = 1,npeak
      if ((energy_c(i) .gt. borneinf) .and. (energy_c(i) .lt. bornesup)) then
         do j = 1,npoint
            absorb(j) = absorb(j) + rot_c(i)/(c**2 + (lambda_c(j)-energy_c(i))**2)
         enddo
      endif
   enddo
   absorb = absorb*q*lambda_c*c

endif

deallocate(energy_c,rot_c,freq_c)

end
