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


write(6,'(A18, i3)') 'Number of state : ', Maxstate
write(6,'(A1)')

! Compute distance between residues
call distance

! Set convolution parameters
call prep_convolution


!--Coupling Zone--!

!Strongest transition
if (stronghest) then
   call cpl_stronghest
endif

!All combination
if (allcomb) then
   call cpl_allcomb
endif

!Extend
if (extend) then
   call cpl_extend
endif

!Superposition without coupling
if (nocoupling) then
   call cpl_nocoupling
endif

deallocate(energy,ground_energy,R,Rn,el_dip)
deallocate(rotatory,freq,lambda_c)


call deal

write(6,'(A1)') " "
write(6,'(A14)') "End of Frenkel"

end
