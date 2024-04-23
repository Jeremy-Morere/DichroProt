subroutine prep_convolution

!Initialize parametres for convolution

use declare
implicit none

!The spectroscopic window and the step are precises in the input file
npoint = floor((bornesup - borneinf)/dw)
dw = (bornesup-borneinf)/npoint

allocate(lambda_c(npoint))
lambda_c = 0.0d0

lambda_c(1) = borneinf
do i = 2, npoint
   lambda_c(i) = lambda_c(i-1) + dw
enddo

!Enlarge the window to take account the foot of near peaks.
bornesup = bornesup + 3*fwhm
borneinf = min(borneinf - 3*fwhm, 0.0d0)

end
