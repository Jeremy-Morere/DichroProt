subroutine store_ss

!The cysteines have to be taken by pair

use declare
implicit none
real*8 :: cnorm

!STEP 1: Find and Store cystine by pair in a temporary array
allocate(res_num_ss(Nss))
allocate(massecenter_ss(Nss,3))
allocate(list_atom_ss(Nss))
res_num_ss = 0
massecenter_ss = 0.0d0
k = 0

do i=1,Nss
 do j=1,i-1
   if ( cnorm(atom_coord(Narom+i,4,:) - atom_coord(Narom+j,4,:)) .le. 3.0 ) then

      res_num_ss(k+1) = res_num(Narom+i)
      massecenter_ss(k+1,:) = massecenter(Narom+i,:)
      list_atom_ss(k+1) = list_atom(Narom+i)

      res_num_ss(k+2) = res_num(Narom+j)
      massecenter_ss(k+2,:) = massecenter(Narom+j,:)
      list_atom_ss(k+2) = list_atom(Narom+j)

      k=k+2
   endif
 enddo
enddo

!STEP 2: Send back in the good order
do i=1,Nss
   res_num(Narom+i) = res_num_ss(i)
   massecenter(Narom+i,:) = massecenter_ss(i,:)
   list_atom(Narom+i) = list_atom_ss(i)
enddo

end
