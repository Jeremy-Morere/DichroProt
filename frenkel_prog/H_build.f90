subroutine H_build

!Build the effectif hamiltonian according to the matrix method.

use declare
implicit none
real*8 :: dip_dip

allocate(H(Nband,Nband))
H = 0.0d0

!Consider only transitions in the near UV.
do t1=1,Nband
 do t2=1,t1
  
  r1=transi_band(t1,1)
  s1=transi_band(t1,2)
  r2=transi_band(t2,1)
  s2=transi_band(t2,2)

  if ((r1 .eq. r2) .and. (s1 .eq. s2)) then !diagonal element
   H(t1,t1) = energy(r1,s1)
  
  elseif (r1 .ne. r2) then !off-diagonal element
   H(t1,t2) = dip_dip(r1,0,s1,r2,0,s2) !Use the dipole-dipole approxiation
   H(t2,t1) = H(t1,t2)
  endif
  
 enddo
enddo

end
