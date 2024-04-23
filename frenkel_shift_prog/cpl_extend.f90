subroutine cpl_extend

!Every transitions of every residus are coupled.

use declare
implicit none

!write(6,'(A37)') "#-Coupling all transitions (extended)"

!Initialization
output = 'extend'
allocate(absorb(npoint))
absorb = 0.0d0

!Coupled part
call H_build
call H_diag(Nresid*Maxstate)
call strenght
call convolution

!!!!!!!!!!!!!!!!!!!
do i=1,Nresid*Maxstate
   if ((energy_c(i) .gt. borneinf) .and. (energy_c(i) .lt. bornesup)) then
      write(6,'(a2,2f10.4)') '| ', 1239.8/energy_c(i), rot_c(i)
   endif
enddo
!!!!!!!!!!!!!!!!!!!

!Writing output
call out_convolution

deallocate(absorb)

end
