subroutine rpdb

!Read the .pdb file.

use declare
implicit none

allocate(atom_coord(Nresid,Maxatom,3))
allocate(atom_charge(Nresid,Maxatom))
atom_coord = 0.0d0
atom_charge = 0.0d0

open(10,file=pdb_file,form='formatted')
do
   read(10,'(A128)',iostat=error) line
   if (error .ne. 0) exit

   if (line(1:4) .eq. "ATOM") then
      read(line(21:26),'(i6)') at
      
      if (Narom .ne. 0) then
      do i=1,Narom  !Aromatic
         if (at .eq. res_num(i)) then
            read(10,'(A128)') line
            read(10,'(A128)') line !Ignor backbone
            read(10,'(A128)') line
            
            do a=1,Natom(i)
               read(10,'(30x,3f8.3,23x,a1)') atom_coord(i,a,1), atom_coord(i,a,2), atom_coord(i,a,3),at_type
               if (at_type .eq. "H") atom_charge(i,a) = 1.0d0 
               if (at_type .eq. "C") atom_charge(i,a) = 6.0d0
               if (at_type .eq. "N") atom_charge(i,a) = 7.0d0
               if (at_type .eq. "O") atom_charge(i,a) = 8.0d0
               if (at_type .eq. "S") atom_charge(i,a) = 1.6d1

            enddo
            read(10,'(A128)') line !Ignor backbone
            read(10,'(A128)') line
         endif
      enddo
      endif

      if (Nss .ne. 0) then
      do i=1, Nss  !Disulfite
      do j=0,1
         k = Narom+i
         if (at .eq. res_num(Narom+2*i+j-1)) then
            read(10,'(A128)') line
            read(10,'(A128)') line !Ignor backbone
            read(10,'(A128)') line

            n_start = 1+4*j
            n_end = 4+4*j

            do a=n_start,n_end
               read(10,'(30x,3f8.3,23x,a1)') atom_coord(k,a,1),atom_coord(k,a,2), atom_coord(k,a,3),at_type
               if (at_type .eq. "H") atom_charge(k,a) = 1.0d0
               if (at_type .eq. "C") atom_charge(k,a) = 6.0d0
               if (at_type .eq. "N") atom_charge(k,a) = 7.0d0
               if (at_type .eq. "O") atom_charge(k,a) = 8.0d0
               if (at_type .eq. "S") atom_charge(k,a) = 1.6d1

            enddo
            read(10,'(A128)') line !Ignor backbone
            read(10,'(A128)') line
         endif
      enddo
      enddo
      endif
   endif
enddo

close(10)

!coversion from A to au
atom_coord = atom_coord*1.8897259886

end
