subroutine cpl_extend

!Every transitions of every residus are coupled.

use declare
implicit none

write(6,'(A37)') "#-Coupling all transitions (extended)"

!Initialization
output = 'extend'
allocate(absorb(npoint))
absorb = 0.0d0

!Coupled part
call H_build
call H_diag(Nresid*Maxstate)
call strenght
call convolution

!Writing output
call out_convolution

deallocate(absorb)

end
