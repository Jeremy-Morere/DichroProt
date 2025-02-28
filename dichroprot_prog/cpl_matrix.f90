subroutine cpl_matrix

!Couple every transition and residue following the matrix method.

use declare
implicit none

write(6,'(A33)') "#-Coupling with the matrix method"
output="matrix"


!FIRST PART: Build and diagonalize the effective Hamiltonian.
call H_build
call H_diag(Nband)

!SECOND PART: Compute rotational strength.
allocate(energy_c(Nband))
allocate(rot_c(Nband))

energy_c = eigval*2.7211d1 !au to eV (2.7211407953d1)
rot_c = 0.0d0

write(6,'(a18,a8,2a12,a8,a12,a8)') '|                 ','Energy', 'Rot', 'mag-el','%','el-el','%'

do ti=1,Nband

 !Initialization of contribution
 mat_mag_el = 0.0d0
 mat_el_el  = 0.0d0
 ri=transi_band(ti,1)
 si=transi_band(ti,2)

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

   mat_mag_el = mat_mag_el - DOT_PRODUCT(mag_dip(r1,s1,:),el_dip(r2,s2,:))*H(t1,ti)*H(t2,ti)*energy(r2,s2)/eigval(ti)
   mat_el_el  = mat_el_el  + DOT_PRODUCT(R(r1,r2,:),cross(:))*H(t1,ti)*H(t2,ti)*energy(r1,s1)*energy(r2,s2)/(4.0d0*eigval(ti))
!   mat_el_el  = mat_el_el  + DOT_PRODUCT(R(r1,r2,:),cross(:))*H(t1,ti)*H(t2,ti)*energy(r1,s1)*energy(r2,s2)/(137*4.0d0*eigval(ti))

  enddo
 enddo

 !conversion from au to 10**-40 esu.erg.cm/G
 mat_mag_el = mat_mag_el*4.714436d2
 mat_el_el  = mat_el_el*4.714436d2

 rot_c(ti) = mat_el_el + mat_mag_el !cgs unit

 !Print rotatory strength and each contribution
 mat = abs(mat_mag_el) + abs(mat_el_el)

 if (ri .le. Narom) then
  write(6,'(a18,f8.2,2e12.3,f8.2,e12.3,f8.2)') &
          '|                 ', 1.239d3/energy_c(ti), rot_c(ti), &
          mat_mag_el, 100*abs(mat_mag_el)/abs(mat), &
          mat_el_el, 100*abs(mat_el_el)/mat
 else
  write(6,'(a18,f8.2,2e12.3,f8.2,e12.3,f8.2)') &
          '|                 ', 1.239d3/energy_c(ti), rot_c(ti), &
          mat_mag_el, 100*abs(mat_mag_el)/abs(mat), &
          mat_el_el, 100*abs(mat_el_el)/mat
 endif
enddo


call convolution

allocate(rot_mat(Nband))
rot_mat = 0.0d0
rot_mat = rot_c


deallocate(H,eigval,energy_c,rot_c)

end
