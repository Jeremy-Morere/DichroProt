subroutine H_build

!Build the Frenkel hamiltonian according to the method chosed.

use declare

!Stronghest and all combinaison
if ((output .eq. "stronghest") .or. (output .eq. "allcomb")) then
allocate(H(Nresid,Nresid))
H = 0.0d0
do r1=1,Nresid
   do r2=1,r1
      if (r1 .eq. r2) then
         H(r1,r1) = energy(r1,comb(r1))
      elseif (Rn(r1,r2) .gt. 1.0d-04) then
         PS1 = DOT_PRODUCT(el_dip(r1,comb(r1),:),el_dip(r2,comb(r2),:))
         PS2 = DOT_PRODUCT(el_dip(r1,comb(r1),:),R(r1,r2,:))
         PS3 = DOT_PRODUCT(el_dip(r2,comb(r2),:),R(r1,r2,:))
         H(r1,r2) = PS1/(Rn(r1,r2)**3) - 3*(PS2 * PS3)/(Rn(r1,r2)**5)
         H(r2,r1) = H(r1,r2)
      endif
   enddo
enddo
endif

!Extend
if (output .eq. "extend") then
allocate(H(Nresid*Maxstate,Nresid*Maxstate))
H = 0.0d0
do r1=1,Nresid
   do r2=1,Nresid!r1
      do s1=1,Maxstate
         do s2=1,Maxstate
   if ((r1 .eq. r2) .and. (s1 .eq. s2)) then !diagonal element
      d = (r1-1)*Maxstate+s1
      H(d,d) = energy(r1,s1)

   elseif (r1 .ne. r2) then !off-diagonal element
      PS1 = DOT_PRODUCT(el_dip(r1,s1,:),el_dip(r2,s2,:))
      PS2 = DOT_PRODUCT(R(r1,r2,:),el_dip(r1,s1,:))
      PS3 = DOT_PRODUCT(R(r1,r2,:),el_dip(r2,s2,:))
      od1 = (r1-1)*Maxstate+s1
      od2 = (r2-1)*Maxstate+s2
      H(od1,od2) = PS1/(Rn(r1,r2)**3) - 3*(PS2 * PS3)/(Rn(r1,r2)**5)
!      H(od2,od1) = H(od1,od2)
   endif
         enddo
      enddo
   enddo
enddo
endif

end
