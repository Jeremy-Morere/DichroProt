program main_geometry

!Get the coordinate of each photo-active side chain

use declare
implicit none
real*8 :: cnorm

call getarg(1,pdb_file)

!STEP 1: Find number of Photo-Activ Residu
call find_number_PAR

!STEP 2: Get atoms coordinate and number
call rpdb

!STEP 3: Compute centre of molecule
call com

!STEP 4: Compute centre of molecule
call store_ss

!STEP 7: Select residu to put in the same qmmm partition

if (.false.) then

write(*,*) "#--Aromatic--#"
do i=1,Narom
   write(*,*) res_num(i)
enddo

write(*,*) "#--SS Bridge--#"

!SS bridge
do i=1,Nss 
 do j=1,i-1
   if ( cnorm(atom_coord(Narom+i,4,:) - atom_coord(Narom+j,4,:)) .le. 3.0 ) then 
      write(*,*) res_num(Narom+i), res_num(Narom+j)
   endif
 enddo
enddo

write(*,*) "#----#"

endif

call woutput


if (.false.) then
!Cette partie permet de detecter les interaction de type pi-stacking


!STEP 4: Compute normal vector of each aromatic 
call normal_vect

!STEP 5: Compute distance between centre
call distance

!STEP 6: Compute angle between com vector and normal vecteur
call compute_angle

!Aromatic 
call min_distance


!Selecton of residus
!Pi stacking 
do i=1,Narom
   do j=1,i-1

      if (Rn(i,j) .le. 3.8) then
         if (angle(i,j,1) .ge. 0.93) then
            if (angle(i,j,2) .ge. 0.93) then
               write(*,*) "face to face"
            endif
         endif
      endif
      if ((Rn(i,j) .le. 5.1) .and. (Rn(i,j) .ge. 5.0))then
         if ((angle(i,j,1) .ge. 0.93) .or. (angle(i,j,1) .le. 0.18)) then
            if (angle(i,j,2) .le. 0.17) then
               write(*,*) "T"
            endif
         endif
      endif
      if ((Rn(i,j) .le. 4.5) .and. (Rn(i,j) .ge. 3.8)) then
         if ((angle(i,j,1) .ge. 0.82) .and. (angle(i,j,1) .le. 0.97)) then
            if (angle(i,j,2) .ge. 0.94) then
               write(*,*) "offset"
            endif
         endif
      endif
   enddo
enddo

endif

call deal

end
