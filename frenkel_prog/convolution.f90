subroutine convolution

!From the list of rotational strengths and transition energies,
!convolute a Lorentzian function.
!The spectroscopic window and the step are specified in the input file.

use declare
implicit none


!Initialization
allocate(absorb(npoint))
absorb = 0.0d0

!---Convolute spectrum---!
q = 1.0d0/(2.294e1*pi) !cf Scott 2021

rot_c = rot_c*q

!Convolute with Lorentzian function
do ti = 1,Nband 
   if ((energy_c(ti) .gt. borneinf) .and. (energy_c(ti) .lt. bornesup)) then
      do j = 1,npoint
         absorb(j) = absorb(j) + rot_c(ti)*lambda_c(j)*fwhm/(fwhm**2 + (lambda_c(j)-energy_c(ti))**2)
      enddo
   endif
enddo

!---Write output---!

write(6,'(a1)')  '|'
write(6,'(a37)') '|Writing the spectrum in a .conv file'
write(6,'(a46)') '|Writing the rotatory strength in a .stck file'
write(6,'(a25)') '#-------------------------------------------'


!Write spectrum
open(14,file='spectrum_'//trim(output)//'.conv',form='formatted',status='replace')
write(14,'(a10,a50)') "E[nm]", "R[L.mol-1.cm-1]"

do i = 1,npoint
   write(14,'(f6.2,f50.5)') 1239.8/lambda_c(i), absorb(i)
enddo

close(14)

!Write rotatory strength
open(15,file='spectrum_'//trim(output)//'.stck',form='formatted',status='replace')
write(15,'(a10,a50)') "E[nm]", "R[L.mol-1.cm-1]"

do ti=1,Nband
   write(15,'(f6.2,f50.5)') 1239.8/energy_c(ti), rot_c(ti)
enddo

close(15)

deallocate(absorb)
end
