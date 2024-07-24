program main_frenkel

use declare
implicit none

!ifort -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O2 -o main_frenkel main_frenkel.f90


write(6,'(a30)')
write(6,'(a30)') '      #----------------------#'
write(6,'(a30)') '      |                      |'
write(6,'(a30)') '      |   Frenkel coupling   |'
write(6,'(a30)') '      |                      |'
write(6,'(a30)') '      #----------------------#'
write(6,'(a30)')  

! Read input file and get gaussian file name
call rinput 
write(6,'(a1)')
write(6,'(a19, i3)') 'Number of residu : ', Nresid
write(6,'(a18, i3)') 'Number of state : ', Nstate
write(6,'(a1)')

! Extract data from gaussian and pdb file
call rlog
call rpdb


! Compute center of charge
call com
! Compute distance between residues
call distance


! Convolution preparation
npoint = floor((bornesup - borneinf)/dw)

write(6,'(a1)') " "
write(6,'(a15,f4.0,a1,f4.0,a3)') "Spectral band: ", borneinf, "-", bornesup, "nm"
tempo = borneinf
borneinf = 1239.8/bornesup !nm to eV
bornesup = 1239.8/tempo    !nm to eV
dw = (bornesup - borneinf)/npoint

write(6,'(a1)') " "


allocate(lambda_c(npoint))
lambda_c = 0.0d0

lambda_c(1) = borneinf
do i = 2, npoint
   lambda_c(i) = lambda_c(i-1) + dw
enddo


!--Coupling Zone--!

!Strongest transition
if (perturbation) then
   call cpl_perturbation
endif

!All combination
if (matrix) then
   call cpl_matrix
endif



deallocate(energy,ground_energy,R,Rn,el_dip)
deallocate(rotatory,freq,lambda_c)


call deal

write(6,'(a1)') " "
write(6,'(a14)') "End of Frenkel"

end
