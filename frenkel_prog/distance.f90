subroutine distance

!Gives center of mass vectors and distance between residus.
!The center is ponderate by atomic number.

use declare
implicit none

!FIRST PART: Compute the center of mass for each residu
call com

!SECOND PART: Compute distance between residu
allocate(Rn(Nresid,Nresid))
allocate(R(Nresid,Nresid,3))
Rn = 0.0d0
R = 0.0d0

do i=1,Nresid
   do j=1,i
!Compute distance vector (Ang to au) 
      R(i,j,1) = (massecenter(j,1)-massecenter(i,1))*1.889725d0
      R(j,i,1) = -R(i,j,1)
      R(i,j,2) = (massecenter(j,2)-massecenter(i,2))*1.889725d0
      R(j,i,2) = -R(i,j,2)
      R(i,j,3) = (massecenter(j,3)-massecenter(i,3))*1.889725d0
      R(j,i,3) = -R(i,j,3)
!Compute norme
      Rn(i,j) = dsqrt(R(i,j,1)**2 + R(i,j,2)**2 + R(i,j,3)**2)
      Rn(j,i) = Rn(i,j)
   enddo
enddo

deallocate(atom_coord,atom_mass)

end
