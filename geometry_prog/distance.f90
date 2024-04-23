subroutine distance

!Compute center of masse vetors of each residus

use declare
implicit none
real*8 :: cnorm

allocate(Rn(Nresid,Nresid))
allocate(R(Nresid,Nresid,3))
Rn = 0.0d0
R = 0.0d0


do i=1,Nresid
   do j=1,i-1
!Compute distance vector
      R(i,j,1) = (massecenter(j,1)-massecenter(i,1))
      R(j,i,1) = -R(i,j,1)
      R(i,j,2) = (massecenter(j,2)-massecenter(i,2))
      R(j,i,2) = -R(i,j,2)
      R(i,j,3) = (massecenter(j,3)-massecenter(i,3))
      R(j,i,3) = -R(i,j,3)
!Compute norme
      Rn(i,j) = cnorm(R(i,j,:))
      Rn(j,i) = Rn(i,j)
   enddo
enddo

end

