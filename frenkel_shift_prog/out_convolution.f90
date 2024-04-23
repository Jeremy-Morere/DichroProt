subroutine out_convolution

!Write spectrum in a .conv file.

use declare
implicit none

write(6,'(a1)')  '|'
write(6,'(a37)') '|Writing the spectrum in a .conv file'
write(6,'(A25)') '#-------------------------------------------'

open(14,file='spectrum_'//trim(output)//'.conv',form='formatted',status='replace')

write(14,'(a10,a50)') "E[nm]", "R[L.mol-1.cm-1]"

do i = 1,npoint
   if (output .eq. "allcomb") then !Ponderation by number of peak plot
      write(14,'(f10.5,f50.5)') 1239.8/lambda_c(i), absorb(i)/Nresid*Maxstate
   else
      write(14,'(f10.5,f50.5)') 1239.8/lambda_c(i), absorb(i)
   endif
enddo

close(14)

end
