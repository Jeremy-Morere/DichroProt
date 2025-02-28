subroutine cpl_perturbation

!Couple every transition and residue following the matrix method.

use declare
implicit none
real*8 :: dip_dip 

write(6,'(A37)') "#-Coupling with the perturbation method"

!Initialization
output="perturbation"
allocate(energy_c(Nband))
allocate(rot_c(Nband))
allocate(Tino_a(Nband))
allocate(Tino_b(Nband))
allocate(Tino_e(Nband))
allocate(Tino_f(Nband))
allocate(Tino(Nband))
energy_c = 0.0d0
rot_c    = 0.0d0
Tino_a = 0.0d0
Tino_b = 0.0d0
Tino_e = 0.0d0
Tino_f = 0.0d0
Tino   = 0.0d0

if (full) then
 write(6,'(a2,3a8,2a12,a8,a12,a8,a12,a8,a12,a8)') '| ','residue','state','Energy','Rot','Tino-a','a %','Tino-b','b %','Tino-e','e %','Tino-f','f %'
else
 write(6,'(a2,3a8,2a12,a8,a12,a8,a12,a8)') '| ','residue','state','Energy','Rot','Tino-a','a %','Tino-b','b %','Tino-f','f %'
endif

do ti=1,Nband
 !Residues and states involved in the transition ti.
 ri=transi_band(ti,1)
 si=transi_band(ti,2)

 !FIRST: Read energy.
 energy_c(ti) = energy(ri,si)*2.7211407953d1 !au to eV

 !SECOND: Compute CD.
 !Tinoco IIIB-22a
 Tino_a(ti) = - DOT_PRODUCT(el_dip(ri,si,:),mag_dip(ri,si,:))

 !Tinoco IIIB-22e
 if (full) then 
  vect(:) = el_dip_perm(ri,si,:) - el_dip_perm(ri,0,:)
  tempo = -DOT_PRODUCT(vect(:),mag_dip(ri,si,:))/energy(ri,si)
  do r1=1,Nresid
   if (r1 .ne. ri) then
    Tino_e(ti) = Tino_e(ti)-tempo*dip_dip(ri,0,si,r1,0,0)
   endif
  enddo
 endif

 do t1=1,Nband
  r1=transi_band(t1,1)
  s1=transi_band(t1,2)
  
  !No coupling between transitions within the same residue or with degenerate transitions.
  if ((r1 .ne. ri) .and. (abs((energy(r1,s1)-energy(ri,si))) .gt. 1.0d-3) ) then

   !Tinoco IIIB-22b
   tempo =-energy(ri,si)*DOT_PRODUCT(el_dip(ri,si,:),mag_dip(r1,s1,:)) &
          -energy(r1,s1)*DOT_PRODUCT(el_dip(r1,s1,:),mag_dip(ri,si,:))
   Tino_b(ti) = Tino_b(ti) &
              -2*tempo*dip_dip(ri,0,si,r1,0,s1)/(energy(r1,s1)**2-energy(ri,si)**2)

   !Tinoco IIIB-22f
   cross(1) = el_dip(r1,s1,2)*el_dip(ri,si,3)-el_dip(r1,s1,3)*el_dip(ri,si,2) 
   cross(2) = el_dip(r1,s1,3)*el_dip(ri,si,1)-el_dip(r1,s1,1)*el_dip(ri,si,3)
   cross(3) = el_dip(r1,s1,1)*el_dip(ri,si,2)-el_dip(r1,s1,2)*el_dip(ri,si,1)
   tempo = DOT_PRODUCT(R(r1,ri,:),cross(:))
   Tino_f(ti) = Tino_f(ti) &
             -dip_dip(ri,0,si,r1,0,s1)*tempo*energy(ri,si)*energy(r1,s1)/(energy(r1,s1)**2-energy(ri,si)**2)
!             -(1.0d0/137)*dip_dip(ri,0,si,r1,0,s1)*tempo*energy(ri,si)*energy(r1,s1)/(energy(r1,s1)**2-energy(ri,si)**2)
  endif

 enddo

enddo
 
 !conversion from au to 10**-40 esu.erg.cm/G
 Tino_a = Tino_a*4.714436d2
 Tino_b = Tino_b*4.714436d2
 Tino_e = Tino_e*4.714436d2
 Tino_f = Tino_f*4.714436d2

 rot_c = Tino_a + Tino_b + Tino_e + Tino_f !cgs unit (10^-40 esu·erg·cm/G)

!Print rotatory strength and each contribution
Tino = abs(Tino_a)+abs(Tino_b)+abs(Tino_e)+abs(Tino_f)

allocate(energy_tempo(Nband))
allocate(rot_pert(Nband)) !Use to compare perturbation and matrix methods in main 
energy_tempo = energy_c
rot_pert     = 0.0d0

do i=1,Nband
 
 ti = MINLOC(energy_tempo,dim=1,mask= energy_tempo > 0)

 energy_tempo(ti) = -1.0d0
 ri=transi_band(ti,1)
 si=transi_band(ti,2)

 !Store rotatory for ratio with matrix
 rot_pert(i) = rot_c(ti)

 if (ri .le. Narom) then
  if (full) then
   write(6,'(a2,2i8,f8.2,2e12.3,f8.2,e12.3,f8.2,e12.3,f8.2,e12.3,f8.2)') &
           '| ', res_num(ri), si, 1.239d3/energy_c(ti), rot_c(ti), &
           Tino_a(ti), 100*abs(Tino_a(ti))/Tino(ti), Tino_b(ti), 100*abs(Tino_b(ti))/Tino(ti), &
           Tino_e(ti), 100*abs(Tino_e(ti))/Tino(ti), Tino_f(ti), 100*abs(Tino_f(ti))/Tino(ti)
  else
   write(6,'(a2,2i8,f8.2,2e12.3,f8.2,e12.3,f8.2,e12.3,f8.2)') &
           '| ', res_num(ri), si, 1.239d3/energy_c(ti), rot_c(ti), &
           Tino_a(ti), 100*abs(Tino_a(ti))/Tino(ti), Tino_b(ti), 100*abs(Tino_b(ti))/Tino(ti), Tino_f(ti), 100*abs(Tino_f(ti))/Tino(ti)
  endif
 
 else
   if (full) then
    write(6,'(a2,2i4,i8,f8.2,2e12.3,f8.2,e12.3,f8.2,e12.3,f8.2,e12.3,f8.2)') &
            '| ', res_num(ri), res_num(ri+1), si, 1.239d3/energy_c(ti), rot_c(ti), &
            Tino_a(ti), 100*abs(Tino_a(ti))/Tino(ti), Tino_b(ti), 100*abs(Tino_b(ti))/Tino(ti), &
            Tino_e(ti), 100*abs(Tino_e(ti))/Tino(ti), Tino_f(ti), 100*abs(Tino_f(ti))/Tino(ti)
  else
   write(6,'(a2,2i4,i8,f8.2,2e12.3,f8.2,e12.3,f8.2,e12.3,f8.2)') &
           '| ', res_num(ri), res_num(ri+1), si, 1.239d3/energy_c(ti), rot_c(ti), &
           Tino_a, 100*abs(Tino_a)/Tino, Tino_b, 100*abs(Tino_b)/Tino, Tino_f, 100*abs(Tino_f)/Tino
  endif
 endif

enddo

call convolution


deallocate(Tino,Tino_a,Tino_b,Tino_e,Tino_f)
deallocate(rot_c,energy_c,energy_tempo)

end
