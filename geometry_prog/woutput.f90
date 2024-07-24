subroutine woutput

!Write an output file use for creating qmmm input

use declare
implicit none


open(11,file="geometry.dat",form='formatted',status='replace')

write(11,'(3i6)') Narom+Nss/2, Narom, Nss/2

do i=1,Narom
   write(11,'(i4,3f16.10,i4,a3,A)') res_num(i), massecenter(i,:), Natom(i),"" ,list_atom(i)
enddo

do i=Narom+1,Nresid,2
   list_atom(i) = TRIM(list_atom(i))//","//TRIM(list_atom(i+1))
   write(11,'(2i4,3f16.10,i4,a3,A)') res_num(i),res_num(i+1),(massecenter(i,:)+massecenter(i+1,:))/2, 4, "", list_atom(i)
enddo

close(11)



end
