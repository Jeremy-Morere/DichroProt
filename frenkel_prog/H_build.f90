subroutine H_build

!Build the Frenkel hamiltonian according to the method chosed.

use declare
implicit none
real*8 :: dip_dip

!Extend
allocate(H(Nband,Nband))
H = 0.0d0

!Take account only the transition in the near UV
do t1=1,Nband
 do t2=1,t1
  
  r1=transi_band(t1,1)
  s1=transi_band(t1,2)
  r2=transi_band(t2,1)
  s2=transi_band(t2,2)

  if ((r1 .eq. r2) .and. (s1 .eq. s2)) then !diagonal element
   H(t1,t1) = energy(r1,s1)
  
  elseif (r1 .ne. r2) then !off-diagonal element
   H(t1,t2) = dip_dip(r1,0,s1,r2,0,s2) !PS1/(Rn(r1,r2)**3) - 3*(PS2 * PS3)/(Rn(r1,r2)**5)
   H(t2,t1) = H(t1,t2)
  endif
  
 enddo
enddo


!~~~~~~~~~~~~~~~~~~~~~~~~
if (.false.) then
open(10,file='H_eff.dat')

write(filename,'(a2,i2,a7)') "'(", Nband, "f10.5)'"

filename=trim(filename)

write(6,*) filename

do t1=1,Nband
   write(10,'(25f15.10)') (H(t1,t2), t2=1,Nband)
enddo

close(10)
endif
!~~~~~~~~~~~~~~~~~~~~~~~~

end
