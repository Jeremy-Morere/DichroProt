subroutine cpl_perturbation

!Couple every transition and residu following the matrix method

use declare
implicit none
real*8 :: dip_dip 

write(6,'(A37)') "#-Coupling with the perturbation method"

!nitialisation
output="perturbation"
allocate(energy_c(Nband))
allocate(rot_c(Nband))

energy_c = 0.0d0
rot_c    = 0.0d0


do ti=1,Nband
 ri=transi_band(ti,1) 
 si=transi_band(ti,2)

 !FIRST : compute energy
 energy_c(ti) = energy(ri,si)*2.7211407953d1 !au to eV

 !SECOND : compute CD
 !Tinoco IIIB-22a
 rot_c(ti) = rot_c(ti) - DOT_PRODUCT(el_dip(ri,si,:),mag_dip(ri,si,:))

 !Tinoco IIIB-22e
 if (full) then 
  vect(:) = el_dip_perm(ri,si,:) - el_dip_perm(ri,0,:)
  tempo = -DOT_PRODUCT(vect(:),mag_dip(ri,si,:))/energy(ri,si)
  do r1=1,Nresid
   if (r1 .ne. ri) then
    rot_c(ti) = rot_c(ti)-tempo*dip_dip(ri,0,si,r1,0,0)
   endif
  enddo
 endif

 do t1=1,Nband
  r1=transi_band(t1,1)
  s1=transi_band(t1,2)

!No coupling between transition of a same residue 
  !if ((r1 .ne. ri) .and. (abs((energy(r1,s1)-energy(ri,si))) .gt. 1.0d-4) ) then
  if (r1 .ne. ri) then

   !Tinoco IIIB-22b
   if (full) then
    tempo =-energy(ri,si)*DOT_PRODUCT(el_dip(ri,si,:),mag_dip(r1,s1,:)) &
           -energy(r1,s1)*DOT_PRODUCT(el_dip(r1,s1,:),mag_dip(ri,si,:))
    rot_c(ti) = rot_c(ti) &
               -2*tempo*dip_dip(ri,0,si,r1,0,s1)/(energy(r1,s1)**2-energy(ri,si)**2)
   endif

   !Tinoco IIIB-22f
   cross(1) = el_dip(r1,s1,2)*el_dip(ri,si,3)-el_dip(r1,s1,3)*el_dip(ri,si,2) 
   cross(2) = el_dip(r1,s1,3)*el_dip(ri,si,1)-el_dip(r1,s1,1)*el_dip(ri,si,3)
   cross(3) = el_dip(r1,s1,1)*el_dip(ri,si,2)-el_dip(r1,s1,2)*el_dip(ri,si,1)
   tempo = DOT_PRODUCT(R(r1,ri,:),cross(:))
   rot_c(ti) = rot_c(ti) &
             -137*dip_dip(ri,0,si,r1,0,s1)*tempo*energy(ri,si)*energy(r1,s1)/(energy(r1,s1)**2-energy(ri,si)**2)

  endif

 enddo
write(6,*) ri, si, rot_c(ti)*4.714436d2
enddo

rot_c = rot_c*4.714436d2 !au to 10**-40 esu.erg.cm/G


!----------------
if (.false.) then
write(6,'(7a15)') 'R', 'Tino a+f', '%','Tino_a', 'Tino_b', 'Tino_e', 'Tino_f'
do ti=1,Nband
 ri=transi_band(ti,1)
 si=transi_band(ti,2)
 Tino_a = 0
 Tino_b = 0
 Tino_e = 0
 Tino_f = 0
 !Tinoco IIIB-22a
 Tino_a = - DOT_PRODUCT(el_dip(ri,si,:),mag_dip(ri,si,:))

 !Tinoco IIIB-22e
 vect(:) = el_dip_perm(ri,si,:) - el_dip_perm(ri,0,:)
 tempo = -DOT_PRODUCT(vect(:),mag_dip(ri,si,:))/energy(ri,si)
 do r1=1,Nresid
  if (r1 .ne. ri) then
   Tino_e = Tino_e -tempo*dip_dip(ri,0,si,r1,0,0)
  endif
 enddo


 do t1=1,Nband
  r1=transi_band(t1,1)
  s1=transi_band(t1,2)

!No coupling with the same residue or if quasi-degenerat
  if ((r1 .ne. ri) .and. (abs((energy(r1,s1)-energy(ri,si))) .gt. 1.0d-4) ) then
   !Tinoco IIIB-22b
   tempo =-energy(ri,si)*DOT_PRODUCT(el_dip(ri,si,:),mag_dip(r1,s1,:)) &
          -energy(r1,s1)*DOT_PRODUCT(el_dip(r1,s1,:),mag_dip(ri,si,:))
   Tino_b = Tino_b -2*tempo*dip_dip(ri,0,si,r1,0,s1)/(energy(r1,s1)**2-energy(ri,si)**2)

   !Tinoco IIIB-22f
   cross(1) = el_dip(r1,s1,2)*el_dip(ri,si,3)-el_dip(r1,s1,3)*el_dip(ri,si,2)
   cross(2) = el_dip(r1,s1,3)*el_dip(ri,si,1)-el_dip(r1,s1,1)*el_dip(ri,si,3)
   cross(3) = el_dip(r1,s1,1)*el_dip(ri,si,2)-el_dip(r1,s1,2)*el_dip(ri,si,1)
   tempo = DOT_PRODUCT(R(r1,ri,:),cross(:))
   Tino_f = Tino_f -137*dip_dip(ri,0,si,r1,0,s1)*tempo*energy(ri,si)*energy(r1,s1)/(energy(r1,s1)**2-energy(ri,si)**2)
  endif

 enddo

write(6,'(7f15.5)') rot_c(ti), (Tino_a+Tino_f)*4.714436d2, 100*(rot_c(ti)-Tino_a-Tino_f)/rot_c(ti), Tino_a*4.714436d2, Tino_b*4.714436d2, Tino_e*4.714436d2, Tino_f*4.714436d2

enddo
endif
!----------------




call convolution

deallocate(rot_c,energy_c)

end
