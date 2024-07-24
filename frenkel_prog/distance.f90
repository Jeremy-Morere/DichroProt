subroutine distance

!Gives center of charge vectors and distance between residus.
!The center is ponderate by atomic number.

use declare
implicit none

allocate(Rn(Nresid,Nresid))
allocate(R(Nresid,Nresid,3))
Rn = 0.0d0
R = 0.0d0

do i=1,Nresid
   do j=1,i
!Compute distance vector 
      R(i,j,1) = (chargecenter(j,1)-chargecenter(i,1))
      R(j,i,1) = -R(i,j,1)
      R(i,j,2) = (chargecenter(j,2)-chargecenter(i,2))
      R(j,i,2) = -R(i,j,2)
      R(i,j,3) = (chargecenter(j,3)-chargecenter(i,3))
      R(j,i,3) = -R(i,j,3)
!Compute norme
      Rn(i,j) = dsqrt(R(i,j,1)**2 + R(i,j,2)**2 + R(i,j,3)**2)
      Rn(j,i) = Rn(i,j)
   enddo
enddo

deallocate(chargecenter)

end
