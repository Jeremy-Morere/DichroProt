!Declaration module
module declare
real*8,allocatable :: absorb(:), lambda(:)
real*8             :: abs_tempo
character*64       :: input,ext,output
integer            :: i,j,N_line,N_frame,error
character*128      :: line
character*1        :: matrix_t, perturbation_t,superposition_t
logical            :: matrix, perturbation,superposition
end module declare



program average_spectrum
!Compute the average spectrum along a trajectory

use declare
implicit none

matrix        = .false.
perturbation  = .false.
superposition = .false.

call getarg(1,input)

!read requested tasks
open(10,file=input)
 do
  read(10,'(A128)',iostat=error) line
  if (error .ne. 0) exit

 if (line(1:14) .eq. 'Perturbation :') then
  read(line(15:128),*) perturbation_t
  if (perturbation_t .eq. 'y') perturbation = .true.
 endif

 if (line(1:8) .eq. 'Matrix :') then
  read(line(9:128),*) matrix_t
  if (matrix_t .eq. 'y') matrix = .true.
 endif

 if (line(1:15) .eq. 'Superposition :') then
  read(line(16:128),*) superposition_t
  if (matrix_t .eq. 'y') superposition = .true.
 endif


 !Read number of frame
 if (line(1:17) .eq. 'Number of frame :') then
  read(line(18:128),*) N_frame
 endif

enddo
close(10)


!Count number of line in the spectrum file
N_line = 0
if (perturbation) then
 open(11,file="perturbation_1.conv")
 read(11,'(A128)',iostat=error) line
 do
  read(11,'(A128)',iostat=error) line
  if (error .ne. 0) exit
  N_line = N_line + 1
 enddo
 close(11)
else
 open(11,file="matrix_1.conv")
 read(11,'(A128)',iostat=error) line
 do
  read(11,'(A128)',iostat=error) line
  if (error .ne. 0) exit
  N_line = N_line + 1
 enddo
 close(11)
endif 


!Compute average 
!If perturbation asked
if (perturbation) then
 output="perturbation"
 call compute_average
endif

!If matrix asked
if (matrix) then
 output="matrix"
 call compute_average
endif

!If superposition asked
if (superposition) then
 output="superposition"
 call compute_average
endif

end


! Current routine to compute the average
subroutine compute_average
use declare
implicit none
allocate(absorb(N_line))
allocate(lambda(N_line))

absorb=0.0d0
lambda=0.0d0
do i=1,N_frame

 !generate file extension
 if (i .lt. 10) then
  write(ext,'(a1,i1,a5)') "_", i, ".conv"
 elseif (i .lt. 100) then
  write(ext,'(a1,i2,a5)') "_", i, ".conv"
  else
  write(ext,'(a1,i3,a5)') "_", i, ".conv"
 endif

 !Read each spectrum and sum absorbtion
 open(12,file=trim(output)//trim(ext))
 read(12,'(A128)') line
 do j=1,N_line
  read(12,'(f6.2,f50.5)') lambda(j), abs_tempo
  absorb(j) = absorb(j) + abs_tempo
 enddo
 close(12)

 !Write output
 open(14,file=trim(output)//'.conv',form='formatted',status='replace')
 do j=1,N_line
  write(14,'(f6.2,f50.5)') lambda(j), absorb(j)/sqrt(N_frame*1.0)
 enddo
 close(14)
enddo

deallocate(lambda,absorb)

end subroutine compute_average
