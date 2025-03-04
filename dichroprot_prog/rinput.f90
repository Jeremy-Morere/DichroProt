subroutine rinput

!Read the input file.

use declare
implicit none

!Initialization
perturbation  = .false.
matrix        = .false.
superposition = .false.
full          = .false. !If .true. compute 

call getarg(1,arg1)


!Read the number of log files from geometry.dat.
open(11,file="geometry.dat")
read(11,'(3i6)') Nresid, Narom, Nss
close(11)

allocate(res_num(Narom+2*Nss))
allocate(namegaus(NResid))

!Read all parameters in the input file
open(10,file=arg1)

do
   read(10,'(A128)',iostat=error) line
   if (error .ne. 0) exit

   !Read log file name and residu number
   if (line(1:14) .eq. '#--Aromatic--#') then
      write(6,'(A18)') 'Residu file name :'
      do i=1,Narom !Aromatic residues
         read(10,*) namegaus(i), res_num(i)
         write(6,*) namegaus(i)
         namegaus(i) = "G09_output/"//trim(namegaus(i))
      enddo
      read(10,'(A128)',iostat=error) line
      do i=0,Nss-1 !SS bridge
         read(10,*) namegaus(Narom+i+1),res_num(Narom+i*2+1),res_num(Narom+i*2+2)
         write(6,*) namegaus(Narom+i+1)
         namegaus(Narom+i+1)="G09_output/"//trim(namegaus(Narom+i+1))
      enddo
   endif

   !Read pdb file name
   if (line(1:5) .eq. 'PDB :') then
      read(line(6:128),*) pdb_file
   endif

   !Read number of state 
   if (line(1:17) .eq. 'Number of state :') then
      read(line(18:128),*) Nstate 
   endif

   !Read if full contribution or not
   if (line(1:19) .eq. 'Full contribution :') then
      read(line(20:128),*) full_t
      if (full_t .eq. 'y') full = .true.
   endif

   !Read requested tasks
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
      if (superposition_t .eq. 'y') superposition = .true.
   endif

 !Read spectra parameter
   !Read spectral window
   if (line(14:28) .eq. 'spectral window') then
      read(10,'(A128)') line
      read(line(14:128),*) borneinf
      read(10,'(A128)') line
      read(line(14:128),*) bornesup
   endif

   !Read spectrum resolution
   if (line(1:12) .eq. 'resolution :') then
      read(line(13:128),*) dw
   endif

   !Read full width at half maximum
   if (line(1:6) .eq. 'fwhm :') then
      read(line,'(6x,f8.2)') fwhm
   endif
   
enddo

close(10)

end
