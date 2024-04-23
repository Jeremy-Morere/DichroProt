subroutine min_distance

!For two residus, compute distance atom per atom and return the minimum

use declare
implicit none
real*8 :: cnorm

allocate(min_distance_list(Narom,Narom))


!Initialisation
do i=1,Narom
 do j=1,i-1
  min_distance_list(i,j) = cnorm((atom_coord(i,1,:) - atom_coord(j,1,:)))
  min_distance_list(j,i) = min_distance_list(i,j) 
 enddo
enddo



do i=1,Narom
 do j=1,Narom-1
  do k=1,Maxat
   do l=1,Maxat
    if ((atom_coord(i,k,1) .ne. 0.0) .or. (atom_coord(j,l,1) .ne. 0.0)) then
     tempo = cnorm((atom_coord(i,k,:) - atom_coord(j,l,:)))

     if (tempo .le. min_distance_list(i,j)) then
      min_distance_list(i,j) = tempo
      min_distance_list(j,i) = min_distance_list(i,j)
     endif
    endif
   enddo
  enddo
 enddo
enddo


end
