subroutine rpdb

!Read pdb and get information about residus and atoms.

use declare
implicit none

!Initialisation residue number
i = 0     !For aromatic
j = Narom !For SS

allocate(atom_coord(Nresid,Maxat,3))
allocate(atom_charge(Nresid,Maxat))
allocate(at_arom_1(Nresid,3)) !use to spot pi-stacking
allocate(at_arom_2(Nresid,3)) !use to spot pi-stacking
allocate(Natom(Nresid))
allocate(res_num(Nresid))
allocate(list_atom(Nresid))

atom_coord = 0.0d0
atom_charge = 0.0d0
at_arom_1 = 0.0d0
at_arom_2 = 0.0d0
Natom = 0
res_num = 0
list_atom = ""

!STEP 2: Get atoms coordinate and number
open(10,file=pdb_file,form='formatted')

do
   read(10,'(A128)',iostat=error) line
   if (error .ne. 0) exit

   Nat = 0

   if (line(18:20) .eq. "PHE") Nat = 14
   if (line(18:20) .eq. "TYR") Nat = 15
   if (line(18:20) .eq. "TRP") Nat = 18
   if (line(18:20) .eq. "CYX") Nat = 4

   !For aromatic
   if (Nat .gt. 10) then 
      i = i + 1
      Natom(i) = Nat

      read(10,'(A128)') line
      read(10,'(A128)') line
      !read residue numbering
      read(10,'(22x,i4)') res_num(i)

      do a=1,Nat
         !read atom numbering, name, coordinate and type
         read(10,'(4x,a7,1x,a4,14x,3f8.3,23x,a1)') at_num, at_name, atom_coord(i,a,:), at_type

         !add atom to the atom list of the residue
         list_atom(i) = TRIM(list_atom(i))//","//TRIM(at_num)

         !Attribute the corresponding atomic number to each atom
         if (at_type .eq. "H") atom_charge(i,a) = 1.0d0
         if (at_type .eq. "C") atom_charge(i,a) = 6.0d0
         if (at_type .eq. "O") atom_charge(i,a) = 8.0d0
         if (at_type .eq. "N") atom_charge(i,a) = 7.0d0
         
         !Get atom coordinates for computing the 
         if (at_name .eq. "CD2") at_arom_1(i,:) = atom_coord(i,a,:)
         if (at_name .eq. "CE2") at_arom_2(i,:) = atom_coord(i,a,:)
      enddo

      read(10,'(A128)') line
      read(10,'(A128)') line     

   endif

   !For ss bridge
   if (Nat .eq. 4) then
      j = j + 1
      Natom(j) = Nat

      read(10,'(A128)') line
      read(10,'(A128)') line
      !read residue numbering
      read(10,'(22x,i4)') res_num(j)

      do a=1,4
         !read atom numbering, coordinate and type
         read(10,'(4x,a7,19x,3f8.3,23x,a1)') at_num, atom_coord(j,a,:), at_type

         !add atom to the atom list of the residue
         list_atom(j) = TRIM(list_atom(j))//","//TRIM(at_num)
         
         !Attribute the corresponding atomic number to each atom
         if (at_type .eq. "H") atom_charge(j,a) = 1.0d0
         if (at_type .eq. "C") atom_charge(j,a) = 6.0d0
         if (at_type .eq. "S") atom_charge(j,a) = 1.6d1
      enddo

      read(10,'(A128)') line
      read(10,'(A128)') line

    endif
enddo

close(10)

!STEP 3: Clean list_atom by removing of blank space

do i=1,Nresid
   call strip_space(list_atom(i))
   list_atom(i) = list_atom(i)(2:)
enddo

end
