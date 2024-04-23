subroutine cpl_allcomb

!Each residu is involved by only one transition.
!Every combinaison of transitions are computed
!and superposed.

use declare
implicit none

!write(6,'(A46)')   "#-Coupling all possible transition combinaison"
!write(6,'(A1)')    "|"
!write(6,'(A30,i)') "|Total number of combinaison: ", (Limstate**Nresid)

!Initialization
output = 'allcomb'
allocate(absorb(npoint))
allocate(comb(Nresid))
absorb = 0.0d0
comb = 1

do 
 if (comb(Nresid) .eq. Limstate+1) exit !Test if is it the last combinaison

  call H_build
  call H_diag(Nresid)
  call strenght
!  call convolution

!Generate next combinaison
  comb(1) = comb(1) + 1
  do i=1,Nresid-1
   if (comb(i) .eq. Limstate+1) then
    comb(i+1) = comb(i+1) + 1
    comb(i) = 1
   endif
  enddo
enddo

!call out_convolution

deallocate(absorb,comb)

end
