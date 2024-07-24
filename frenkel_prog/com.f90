subroutine com

!Compute center of each residus.
!The center is ponderate by atomic number.

use declare
implicit none

!Initialization
allocate(chargecenter(Nresid,3))
chargecenter = 0.0d0

do i=1,Nresid !For each rediues
 !Compute the sum of atomic number
   Charge = 0.0d0
   do a=1,Natom(i)
      Charge = Charge + atom_charge(i,a)

   enddo
 !Compute center without ponderation
   do a=1,Natom(i)
      chargecenter(i,1) = chargecenter(i,1) + atom_charge(i,a)*atom_coord(i,a,1)
      chargecenter(i,2) = chargecenter(i,2) + atom_charge(i,a)*atom_coord(i,a,2)
      chargecenter(i,3) = chargecenter(i,3) + atom_charge(i,a)*atom_coord(i,a,3)
   enddo
 !Ponderate by the sum
   chargecenter(i,:) = chargecenter(i,:)/Charge
enddo

deallocate(atom_coord,atom_charge)

end
