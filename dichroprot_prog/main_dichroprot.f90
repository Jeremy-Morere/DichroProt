program main_dichroprot

use declare
implicit none


write(6,'(A30)')
write(6,'(A30)') '      #----------------------#'
write(6,'(A30)') '      |                      |'
write(6,'(A30)') '      |      Dichroprot      |'
write(6,'(A30)') '      |                      |'
write(6,'(A30)') '      #----------------------#'
write(6,'(A30)')  

! Read input file and get gaussian file name
call rinput 
write(6,'(A19, i3)') 'Number of residu : ', Nresid

! Extract data from gaussian and pdb file
call rlog
call rpdb


write(6,'(A18, i3)') 'Number of state : ', Nstate
write(6,'(A1)')

! Compute center of masse
call coc
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

!Enlarge the window to take account the foot of near peaks.
bornesup = bornesup + 3*fwhm
borneinf = max(borneinf - 3*fwhm, 0.0d0)

!--Coupling Zone--!
!Uncoupled spectra
if (superposition) then
   call cpl_superposition
endif

!Perturbation method
if (perturbation) then
   call cpl_perturbation
endif

!Matrix method
if (matrix) then
   call cpl_matrix
endif


!--Spectre UV
allocate(energy_c(Nband))
allocate(rot_c(Nband))

output="UV"

do ti=1,Nband
 !Residues and states involved in the transition ti.
 ri=transi_band(ti,1)
 si=transi_band(ti,2)

 energy_c(ti) = energy(ri,si)*2.7211407953d1 !au to eV
 rot_c(ti)    = freq(ri,si)
enddo

call convolution

deallocate(energy_c,rot_c)
!--End Spectre UV





deallocate(energy,ground_energy,R,Rn,el_dip,mag_dip)
deallocate(rotatory,freq,lambda_c)


call deal

write(6,'(A1)') " "
write(6,'(A14)') "End of DichroProt"

end
