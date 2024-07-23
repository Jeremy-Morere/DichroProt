subroutine com

!Compute center of masse od each residus

use declare
implicit none

allocate(massecenter(Nresid,3))
massecenter = 0.0d0

do i=1,Nresid
   charge = 0.0d0
   do a=1,Natom(i)
      charge = charge + atom_charge(i,a)
   enddo
   do a=1,Natom(i)
      massecenter(i,1) = massecenter(i,1) + atom_charge(i,a)*atom_coord(i,a,1)/charge
      massecenter(i,2) = massecenter(i,2) + atom_charge(i,a)*atom_coord(i,a,2)/charge
      massecenter(i,3) = massecenter(i,3) + atom_charge(i,a)*atom_coord(i,a,3)/charge
   enddo
enddo

deallocate(atom_charge)

end
