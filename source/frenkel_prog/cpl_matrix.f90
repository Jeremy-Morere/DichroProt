subroutine cpl_matrix

!Couple every transition and residu following the matrix method

use declare
implicit none

write(6,'(A33)') "#-Coupling with the matrix method"
output="matrix"


!Coupled part
call H_build
call H_diag(Nband)

!Computation of rotational strength
allocate(energy_c(Nband))
allocate(rot_c(Nband))

energy_c = eigval*2.7211d1 !au to eV (2.7211407953d1)
rot_c = 0.0d0


do ti=1,Nband 
 do t1=1,Nband
  do t2=1,Nband

   r1=transi_band(t1,1)
   s1=transi_band(t1,2)
   r2=transi_band(t2,1)
   s2=transi_band(t2,2)

!Cross product of electric transition dipoles
   cross(1) = el_dip(r1,s1,2)*el_dip(r2,s2,3)-el_dip(r1,s1,3)*el_dip(r2,s2,2)
   cross(2) = el_dip(r1,s1,3)*el_dip(r2,s2,1)-el_dip(r1,s1,1)*el_dip(r2,s2,3)
   cross(3) = el_dip(r1,s1,1)*el_dip(r2,s2,2)-el_dip(r1,s1,2)*el_dip(r2,s2,1)

   rot_c(ti) = rot_c(ti) + 137*DOT_PRODUCT(R(r1,r2,:),cross(:))*H(t1,ti)*H(t2,ti)*energy(r1,s1)*energy(r2,s2)
   rot_c(ti) = rot_c(ti) - DOT_PRODUCT(mag_dip(r1,s1,:),el_dip(r2,s2,:))*H(t1,ti)*H(t2,ti)*energy(r2,s2)

  enddo
 enddo
 rot_c(ti) = rot_c(ti)/eigval(ti)
enddo

rot_c = rot_c*4.714436d2 !au to 10**-40 esu.erg.cm/G




!----------------------
if (.false.) then
write(6,'(3a15)') 'R','dip x dip', 'mag x dip'
do ti=1,Nband
Tino_a = 0
Tino_b = 0

 do t1=1,Nband
  do t2=1,Nband

   r1=transi_band(t1,1)
   s1=transi_band(t1,2)
   r2=transi_band(t2,1)
   s2=transi_band(t2,2)

!Cross product of electric transition dipoles
   cross(1) = el_dip(r1,s1,2)*el_dip(r2,s2,3)-el_dip(r1,s1,3)*el_dip(r2,s2,2)
   cross(2) = el_dip(r1,s1,3)*el_dip(r2,s2,1)-el_dip(r1,s1,1)*el_dip(r2,s2,3)
   cross(3) = el_dip(r1,s1,1)*el_dip(r2,s2,2)-el_dip(r1,s1,2)*el_dip(r2,s2,1)

   Tino_a = Tino_a + 137*DOT_PRODUCT(R(r1,r2,:),cross(:))*H(t1,ti)*H(t2,ti)*energy(r1,s1)*energy(r2,s2)
   Tino_b = Tino_b - DOT_PRODUCT(mag_dip(r1,s1,:),el_dip(r2,s2,:))*H(t1,ti)*H(t2,ti)*energy(r2,s2)

  enddo
 enddo
write(6,'(3f15.5)') rot_c(ti), Tino_a*4.714436d2, Tino_b*4.714436d2
enddo

endif

!----------------------




call convolution

deallocate(H,eigval)

end
