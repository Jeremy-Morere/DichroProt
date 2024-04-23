subroutine strenght

!Compute rotatory and oscilator strenghts of the coupled system.

use declare
implicit none

if ((output .eq. "stronghest") .or. (output .eq. "allcomb")) then

allocate(energy_c(Nresid))
allocate(rot_c(Nresid))
allocate(freq_c(Nresid))
energy_c = eigval*2.7211407953d1 !au to eV
rot_c = 0.0d0
freq_c = 0.0d0

do i=1,Nresid
   do r1 = 1,Nresid
    do r2 = 1,Nresid
       s1 = comb(s1)
       s2 = comb(r2)

       cross(1) = el_dip(r2,s2,2)*el_dip(r1,s1,3) - el_dip(r2,s2,3)*el_dip(r1,s1,2)
       cross(2) = el_dip(r2,s2,3)*el_dip(r1,s1,1) - el_dip(r2,s2,1)*el_dip(r1,s1,3)
       cross(3) = el_dip(r2,s2,1)*el_dip(r1,s1,2) - el_dip(r2,s2,2)*el_dip(r1,s1,1)

       vect(:) = (massecenter(r2,:)*energy(r2,s2) - massecenter(r1,:)*energy(r1,s1))/(4.556177d5*0.53) !hartree to rb conversion

       freq_c(i) = freq_c(i) + DOT_PRODUCT(el_dip(r1,s1,:),el_dip(r2,s2,:)) *H(r1,i)*H(r2,i)
       rot_c(i) = rot_c(i) + DOT_PRODUCT(vect(:), cross(:)) * H(r1,i)*H(r2,i)

    enddo
   enddo
   rot_c(i) = -rot_c(i)*(eigval(i)/5.481436d2) * 4.714436d2 !au to 10**-40 esu.erg.cm/G 
enddo


else if (output .eq. "extend") then

allocate(energy_c(Nresid*Maxstate))
allocate(rot_c(Nresid*Maxstate))
allocate(freq_c(Nresid*Maxstate))
energy_c = eigval*2.7211407953d1 !au to eV
rot_c = 0.0d0
!freq_c = 0.0d0

do i=1,Nresid*Maxstate
   do r1 = 1,Nresid
    do r2 = 1,Nresid
     do s1 = 1,Maxstate
      do s2 = 1,Maxstate

         cross(1) = el_dip(r2,s2,2)*el_dip(r1,s1,3) - el_dip(r2,s2,3)*el_dip(r1,s1,2)
         cross(2) = el_dip(r2,s2,3)*el_dip(r1,s1,1) - el_dip(r2,s2,1)*el_dip(r1,s1,3)
         cross(3) = el_dip(r2,s2,1)*el_dip(r1,s1,2) - el_dip(r2,s2,2)*el_dip(r1,s1,1)

         vect(:) = (massecenter(r2,:)*energy(r2,s2) - massecenter(r1,:)*energy(r1,s1))/(4.556177d5*0.53) !hartree to rb conversion

         od1 = (r1-1)*Maxstate+s1
         od2 = (r2-1)*Maxstate+s2

         !freq_c(i) = freq_c(i) + DOT_PRODUCT(el_dip(r1,s1,:),el_dip(r2,s2,:)) * H(od1,i)*H(od2,i)
         rot_c(i) = rot_c(i) + DOT_PRODUCT(vect(:), cross(:)) * H(od1,i)*H(od2,i)

      enddo
     enddo
    enddo
   enddo
   rot_c(i) = rot_c(i) * 4.714436d2 !au to 10**-40 esu.erg.cm/G
   write(*,*) rot_c(i), energy_c(i)
enddo

endif

deallocate(H,eigval)

end
