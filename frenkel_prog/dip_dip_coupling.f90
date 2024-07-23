function dip_dip(res1,s1i,s1f,res2,s2i,s2f) result(cpl)

!Compute the coupling between two transition of two different residue
use declare

implicit none
integer, intent (in) :: res1, s1i, s1f, res2, s2i, s2f
real*8               :: PS1, PS2, PS3
real*8               :: cpl

cpl = 0

if (res1 .ne. res2) then !Any coupling in the same residue

   if ((s1i .eq. s1f) .and. (s2i .eq. s2f)) then !1 and 2 with permanent dipole
      PS1 = DOT_PRODUCT(el_dip_perm(res1,s1f,:),el_dip_perm(res2,s2f,:))
      PS2 = DOT_PRODUCT(R(res1,res2,:),el_dip_perm(res1,s1f,:))
      PS3 = DOT_PRODUCT(R(res1,res2,:),el_dip_perm(res2,s2f,:))

   else if (s1i .eq. s1f) then !1 with permanent, 2 with transition
      PS1 = DOT_PRODUCT(el_dip_perm(res1,s1f,:),el_dip(res2,s2f,:))
      PS2 = DOT_PRODUCT(R(res1,res2,:),el_dip_perm(res1,s1f,:))
      PS3 = DOT_PRODUCT(R(res1,res2,:),el_dip(res2,s2f,:))

   else if (s2i .eq. s2f) then !1 with transition, 2 with permanent
      PS1 = DOT_PRODUCT(el_dip(res1,s1f,:),el_dip(res2,s2f,:))
      PS2 = DOT_PRODUCT(R(res1,res2,:),el_dip(res1,s1f,:))
      PS3 = DOT_PRODUCT(R(res1,res2,:),el_dip_perm(res2,s2f,:))

   else !1 with transition, 2 with transition
      PS1 = DOT_PRODUCT(el_dip(res1,s1f,:),el_dip(res2,s2f,:))
      PS2 = DOT_PRODUCT(R(res1,res2,:),el_dip(res1,s1f,:))
      PS3 = DOT_PRODUCT(R(res1,res2,:),el_dip(res2,s2f,:))
   endif

   cpl = PS1/(Rn(res1,res2)**3) - 3*(PS2 * PS3)/(Rn(res1,res2)**5)

endif

end function

