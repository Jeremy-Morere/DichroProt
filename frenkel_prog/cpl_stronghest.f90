subroutine cpl_stronghest

!Keep the stronghest transition and couple them.

use declare
implicit none

write(6,'(A33)') "#-Coupling stronghest transitions"
write(6,'(A1)')  "|"

output = 'stronghest'
allocate(absorb(npoint))
allocate(comb(Nresid))
absorb = 0.0d0
comb = 1

!Find transitions with strongest rotatory power
do i=1,Nresid
   do j=1,Maxstate
      if (abs(rotatory(i,j)) .gt. abs(rotatory(i,comb(i)))) comb(i) = j
   enddo
enddo


write(6,'(A48)')            "|The stronghest transition for each residu are :"
do i=1,Nresid
 write(6,'(A10,i3,A13,i3)') "|  Residu ", i, ": transition ", comb(i)
enddo


call H_build
call H_diag(Nresid)
call strenght

call convolution
call out_convolution

deallocate(absorb,comb)

end
