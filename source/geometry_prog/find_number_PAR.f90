subroutine find_number_PAR

!Find phot-actif residus

use declare
implicit none

!Initialisation
Narom = 0
Nss = 0
Maxat = 0

open(10,file=pdb_file,form='formatted')

!STEP 1: Find number of Photo-Activ Residu

do
   read(10,'(A128)',iostat=error) line
   if (error .ne. 0) exit
!Find Phenylalanine
   if (line(18:20) .eq. "PHE") then
     Narom = Narom + 1
     do i=1,19
        read(10,'(A128)',iostat=error) line
     enddo
     if (Maxat .lt. 14) Maxat = 14
   endif
!Find Tyrosine
   if (line(18:20) .eq. "TYR") then
     Narom = Narom + 1
     do i=1,20
        read(10,'(A128)',iostat=error) line
     enddo
     if (Maxat .lt. 15) Maxat = 15
   endif
!Find Troptophane
   if (line(18:20) .eq. "TRP") then
     Narom = Narom + 1
     do i=1,23
        read(10,'(A128)',iostat=error) line
     enddo
     if (Maxat .lt. 18) Maxat = 18
   endif
!Find cysteine involved in SS-Bridge
   if (line(18:20) .eq. "CYX") then
     Nss = Nss + 1
     do i=1,9
        read(10,'(A128)',iostat=error) line
     enddo
     if (Maxat .lt. 4) Maxat = 4
   endif
enddo


close(10)

Nresid = Narom + Nss



end
