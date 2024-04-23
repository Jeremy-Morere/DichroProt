program main_geometry

!Get the coordinate of each photo-active side chain

use declare
implicit none

!ifort -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O2 -o main_frenkel main_frenkel.f90

!Initialisation
Narom = 0
Nss = 0

call getarg(1,pdb_file)

open(10,file=pdb_file,form='formatted') 

!STEP 1: Find number of Photo-Activ Residu

do
   read(10,'(A128)',iostat=error) line
   if (error .ne. 0) exit
!Find Phenylalanine
   if (line(18:20) .eq. "PHE") then
     Narom = Narom + 1
     do i=1,19
        read(10,'(A128)',iostat=error) line
     enddo
   endif
!Find Tyrosine
   if (line(18:20) .eq. "TYR") then
     Narom = Narom + 1
     do i=1,20
        read(10,'(A128)',iostat=error) line
     enddo
   endif
!Find Troptophane
   if (line(18:20) .eq. "TRP") then
     Narom = Narom + 1
     do i=1,23
        read(10,'(A128)',iostat=error) line
     enddo
   endif
!Find cysteine involved in SS-Bridge
   if (line(18:20) .eq. "CYX") then
     Narom = Narom + 1
     do i=1,9
        read(10,'(A128)',iostat=error) line
     enddo
   endif
enddo

rewind(10)

Nresid = Narom + Nss

!Initialisation residu number
i = 0 !For aromatic
j = 0 !For SS
k = 0 !For every one

allocate(atom_coord(Nresid,18,3))
atom_coord = 0.0d0
allocate(atom_charge(Nresid,18))
atom_charge = 0.0d0
allocate(res_arom_num(Narom))
allocate(res_ss_num(Nss))
res_arom_num = 0
res_ss_num = 0
allocate(at_arom_1(Nresid,3))
allocate(at_arom_2(Nresid,3))
at_arom_1 = 0.0d0
at_arom_2 = 0.0d0
!allocate(list_atom(Nresid))

!STEP 2: Get atoms coordinate and number
do
   read(10,'(A128)',iostat=error) line
   if (error .ne. 0) exit

   Nat = 0

   if (line(18:20) .eq. "PHE") Nat = 14
   if (line(18:20) .eq. "TYR") Nat = 15
   if (line(18:20) .eq. "TRP") Nat = 18
   if (line(18:20) .eq. "CYX") Nat = 4

   if (Nat .gt. 10) then  !For aromatic 
      i = i + 1
      k = k + 1

      read(10,'(A128)') line
      read(10,'(A128)') line
      read(10,'(22x,i4)') res_arom_num(i)

      do a=1,Nat
         read(10,'(4x,a7,a6,13x,3f8.3,23x,a1)') at_num, at_name,atom_coord(k,a,1), atom_coord(k,a,2), atom_coord(k,a,3),at_type

         if (at_type .eq. "H") atom_charge(k,a) = 1.0d0
         if (at_type .eq. "C") atom_charge(k,a) = 6.0d0
         if (at_type .eq. "O") atom_charge(k,a) = 8.0d0
         if (at_type .eq. "N") atom_charge(k,a) = 7.0d0

         if (at_name .eq. "CD2") at_arom_1(i,:) = atom_coord(i,a,:)
         if (at_name .eq. "CE2") at_arom_2(i,:) = atom_coord(i,a,:)
      enddo

      read(10,'(A128)') line
      read(10,'(A128)') line

   endif

   if (Nat .eq. 4) then !For ss bridge
      j = j + 1
      k = k + 1

      read(10,'(A128)') line
      read(10,'(A128)') line
      read(10,'(22x,i4)') res_ss_num(j)

      do a=1,4
         read(10,'(4x,a7,a6,13x,3f8.3,23x,a1)') at_num, at_name,atom_coord(k,a,1), atom_coord(k,a,2), atom_coord(k,a,3),at_type

         if (at_type .eq. "H") atom_charge(k,a) = 1.0d0
         if (at_type .eq. "C") atom_charge(k,a) = 6.0d0
         if (at_type .eq. "S") atom_charge(k,a) = 1.6d1
      enddo

      read(10,'(A128)') line
      read(10,'(A128)') line

    endif
enddo


close(10)

!Concatenate res_arom_num and res_ss_num
allocate(res_num(Nresid))

do i=1,Narom
   res_num(i) = res_num(i)
enddo
do j=1,Nss
   res_num(Narom+j) = res_num(j)
enddo

write(*,*) "STEP 2 finish"

!STEP 3: Compute centre of molecule
allocate(massecenter(Nresid,3))
massecenter = 0.0d0


do i=1,Nresid
   charge = 0.0d0
   do a=1,18
      charge = charge + atom_charge(i,a)
   enddo
   do a=1,18
      massecenter(i,1) = massecenter(i,1) + atom_charge(i,a)*atom_coord(i,a,1)/charge
      massecenter(i,2) = massecenter(i,2) + atom_charge(i,a)*atom_coord(i,a,2)/charge
      massecenter(i,3) = massecenter(i,3) + atom_charge(i,a)*atom_coord(i,a,3)/charge
   enddo
enddo

write(*,*) "STEP 3 finish"

!STEP 4: Compute normal vector of each aromatic 
allocate(normal_vector(3,Narom))

do i=1,Narom
   vect_1(:) = at_arom_1(i,:) - massecenter(i,:)
   vect_2(:) = at_arom_2(i,:) - massecenter(i,:)

   normal_vector(i,1) = vect_1(2)*vect_2(3) - vect_1(3)*vect_2(2)
   normal_vector(i,2) = vect_1(3)*vect_2(1) - vect_1(1)*vect_2(3)
   normal_vector(i,3) = vect_1(1)*vect_2(2) - vect_1(2)*vect_2(1)
   
   norm  = dsqrt(normal_vector(i,1)**2 + normal_vector(i,2)**2 + normal_vector(i,3)**2)
   normal_vector(i,:) = normal_vector(i,:)/norm
enddo
write(*,*) "STEP 4 finish"


!STEP 5: Compute distance between centre
if (allocated(R)) deallocate(R)
if (allocated(Rn)) deallocate(Rn)
write(*,*) "J'ai bien dealloué"
allocate(Rn(Nresid,Nresid))
allocate(R(Nresid,Nresid,3))
write(*,*) "Allocate finish"

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
      Rn(i,j) = dsqrt(R(i,j,1)**2 + R(i,j,2)**2 + R(i,j,3)**2)
      Rn(j,i) = Rn(i,j)

   enddo
enddo

!STEP 6: Compute angle between com vector and normal vecteur
allocate(angle(Narom,Narom,2))

do i=1,Narom
   do j=1,Narom
      angle(i,j,1) = ABS(DOT_PRODUCT(normal_vector(i,:),R(i,j,:))/Rn(i,j))
      angle(i,j,2) = ABS(DOT_PRODUCT(normal_vector(j,:),R(i,j,:))/Rn(i,j))
   enddo
enddo






do i=1,Nresid
write(*,'(i5,6f8.3)') res_num(i), at_arom_1(i,:), at_arom_2(i,:) 
enddo

call deal

open(11,file="distance_table.dat",form='formatted',status='replace')
close(11)

end
