subroutine compute_angle

!Compute Angle between:
!1: normal vector or vectors between center of charge 
!2: normal vectors

use declare
implicit none

allocate(angle(Narom,Narom,2))

do i=1,Narom
   do j=1,i-1
      angle(i,j,1) = ABS(DOT_PRODUCT(normal_vector(i,:),R(i,j,:))/Rn(i,j))
      angle(i,j,2) = ABS(DOT_PRODUCT(normal_vector(j,:),normal_vector(i,:)))
      angle(j,i,1) =  angle(i,j,1)
      angle(j,i,2) =  angle(i,j,2)

   enddo
enddo


end

