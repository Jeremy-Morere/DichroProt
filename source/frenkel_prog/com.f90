subroutine com

!Compute center of each residus.
!The center is ponderate by atomic number.

use declare
implicit none

!Initialization
allocate(massecenter(Nresid,3))
massecenter = 0.0d0

do i=1,Nresid !For each rediues
 !Compute the sum of atomic number
   Mass = 0.0d0
   do a=1,Natom(i)
      Mass = Mass + atom_mass(i,a)

   enddo
 !Compute center without ponderation
   do a=1,Natom(i)
      massecenter(i,1) = massecenter(i,1) + atom_mass(i,a)*atom_coord(i,a,1)
      massecenter(i,2) = massecenter(i,2) + atom_mass(i,a)*atom_coord(i,a,2)
      massecenter(i,3) = massecenter(i,3) + atom_mass(i,a)*atom_coord(i,a,3)
   enddo
 !Ponderate by the sum
   massecenter(i,:) = massecenter(i,:)/Mass
enddo



end
