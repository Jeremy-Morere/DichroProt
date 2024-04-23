subroutine cpl_nocoupling

!Superposition of CD spectra of each residu.

use declare
implicit none

write(6,'(A32)') "#-Superposition without coupling"

!Initialization
output = 'no_coupling'
allocate(absorb(npoint))
absorb = 0.0d0

allocate(energy_c(Nresid*Maxstate))
allocate(rot_c(Nresid*Maxstate))
allocate(freq_c(Nresid*Maxstate))

!Read data in general list to the convolution list
do i=1,Nresid
 do j=1,Maxstate
  energy_c(i+(j-1)*Nresid) = energy(i,j)*2.7211407953d1 !au to eV
  rot_c(i+(j-1)*Nresid) = rotatory(i,j)
  freq_c(i+(j-1)*Nresid) = freq(i,j) 
 enddo
enddo

!Write peak in screen
write(6,'(a22)') '|   E[nm]       R[csg]'
do i=1,Nresid*Maxstate
   if ((energy_c(i) .gt. borneinf) .and. (energy_c(i) .lt. bornesup)) then
      write(6,'(a2,2f10.4)') '| ', 1239.8/energy_c(i), rot_c(i)
   endif
enddo

call convolution
call out_convolution

deallocate(absorb)

end
