subroutine rgeo

!Read geometry.dat file.

use declare
implicit none

!FIRST PART: Read numbers of atom and center of charge for each residue
allocate(Natom(Nresid))
Natom = 0

open(11,file="geometry.dat",form='formatted')
read(11,'(3i6)') Nresid, Narom, Nss


allocate(chargecenter(Nresid,3))
chargecenter = 0.0d0

!Read aromatic center of charge and number of atome
do i=1,Narom
 read(10,'(4x,3f16.10,i4)') chargecenter(i,1), chargecenter(i,2), chargecenter(i,3), Natom(i)
enddo

!Read ss bridge center of charge and number of atome
do i=1,Nss
 read(10,'(8x,3f16.10,i4)') chargecenter(i,1), chargecenter(i,2), chargecenter(i,3), Natom(i)
enddo

close(11)

end
