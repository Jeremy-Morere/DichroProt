subroutine normal_vect

!Compute center of masse od each residus

use declare
implicit none
real*8 :: cnorm

allocate(normal_vector(Narom,3))

do i=1,Narom
   vect_1(:) = at_arom_1(i,:) - massecenter(i,:)
   vect_2(:) = at_arom_2(i,:) - massecenter(i,:)

   normal_vector(i,1) = vect_1(2)*vect_2(3) - vect_1(3)*vect_2(2)
   normal_vector(i,2) = vect_1(3)*vect_2(1) - vect_1(1)*vect_2(3)
   normal_vector(i,3) = vect_1(1)*vect_2(2) - vect_1(2)*vect_2(1)

   norm  = cnorm(normal_vector(i,:))
   normal_vector(i,:) = normal_vector(i,:)/norm
enddo

deallocate(at_arom_1,at_arom_2)

end

