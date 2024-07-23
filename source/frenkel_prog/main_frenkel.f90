program main_frenkel

use declare
implicit none

!ifort -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O2 -o main_frenkel main_frenkel.f90


write(6,'(A30)')
write(6,'(A30)') '      #----------------------#'
write(6,'(A30)') '      |                      |'
write(6,'(A30)') '      |   Frenkel coupling   |'
write(6,'(A30)') '      |                      |'
write(6,'(A30)') '      #----------------------#'
write(6,'(A30)')  

! Read input file and get gaussian file name
call rinput 
write(6,'(A19, i3)') 'Number of residu : ', Nresid

! Extract data from gaussian and pdb file
call rlog
call rpdb


write(6,'(A18, i3)') 'Number of state : ', Nresid
write(6,'(A1)')

! Compute center of masse
call com
! Compute distance between residues
call distance

!--Coupling Zone--!

! Convolution preparation
allocate(lambda_c(npoint))
lambda_c = 0.0d0

lambda_c(1) = borneinf
do i = 2, npoint
   lambda_c(i) = lambda_c(i-1) + dw
enddo


!Strongest transition
if (perturbation) then
   call cpl_perturbation
endif

!All combination
if (matrix) then
   call cpl_matrix
endif

!deallocate(energy,ground_energy,R,Rn,el_dip)
!deallocate(rotatory,freq,lambda_c)


call deal

write(6,'(A1)') " "
write(6,'(A14)') "End of Frenkel"

end
