subroutine rinput

!Read the input file.

use declare
implicit none

stronghest = .false.
allcomb = .false.
extend = .false.
nocoupling = .false.

call getarg(1,arg1)

open(10,file=arg1)

!read input and deduce the number of log file
do
   read(10,'(A128)',iostat=error) line
   if (error .ne. 0) exit
   
   if (line(1:14) .eq. '#--Aromatic--#') then
      Nresid=0
      do
         read(10,'(A128)') line
         if (line(1:6) .eq. '#----#') exit
         if (line(1:15) .eq. "#--SS Bridge--#") Narom = Nresid
         Nresid = Nresid + 1
      enddo
   exit
   endif
enddo
Nresid = Nresid - 1
Nss = Nresid - Narom

allocate(res_num(Narom+2*Nss))
allocate(namegaus(NResid))

rewind 10

do
   read(10,'(A128)',iostat=error) line
   if (error .ne. 0) exit

   !read log file name and residu number
   if (line(1:14) .eq. '#--Aromatic--#') then
      write(6,'(A18)') 'Residu file name :'
      do i=1,Narom
         read(10,*) namegaus(i), res_num(i)
         namegaus(i) = trim(namegaus(i))
         write(6,*) namegaus(i)
      enddo
      read(10,'(A128)',iostat=error) line
      do i=0,Nss-1
         read(10,*) namegaus(Narom+i+1),res_num(Narom+i*2+1),res_num(Narom+i*2+2)
         namegaus(Narom+i+1) = trim(namegaus(Narom+i+1))
         write(6,*) namegaus(Narom+i)
      enddo
   endif
   
   !read pdb file name
   if (line(1:5) .eq. 'PDB :') then
      read(line(6:128),*) pdb_file
   endif

   !read prmtop file name
   if (line(1:8) .eq. 'PRMTOP :') then
      read(line(9:128),*) prmtop_file
   endif


   !read requested tasks
   if (line(1:17) .eq. '#Type of coupling') then
      read(10,'(A128)') line
      if (line(25:25) .eq. 'y') then
         stronghest = .true.
      endif

      read(10,'(A128)') line
      if (line(25:25) .eq. 'y') then
         allcomb = .true.
      endif

      read(10,'(A128)') line
      if (line(25:25) .eq. 'y') then
         extend = .true.
      endif

      read(10,'(A128)') line
      if (line(25:25) .eq. 'y') then
         nocoupling = .true.
      endif
   endif
   
   !If allcomb, read number of excited state to take in calculation
   if ((line(1:22) .eq. "number excited state :") .and. (allcomb)) then
      read(line(23:128),*) Limstate   
   endif

   !read spectral window
   if (line(14:28) .eq. 'spectral window') then
      read(10,'(A128)') line
      read(line(14:128),*) borneinf
      read(10,'(A128)') line
      read(line(14:128),*) bornesup
   endif

   !read spectrum resolution
   if (line(1:12) .eq. 'resolution :') then
      read(line(13:128),*) dw
   endif

   !read full width at half maximum
   if (line(1:6) .eq. 'fwhm :') then
      read(line,'(6x,f8.2)') fwhm
   endif
   
   !read convolution fonction 
   if (line(1:10) .eq. 'function :') then
      read(line,'(11x,a1)') type
   endif

enddo

close(10)


write(6,'(a1)') " "
write(6,'(a15,f4.0,a1,f4.0,a3)') "Spectral band: ", borneinf, "-", bornesup, " nm"
tempo = borneinf
borneinf = 1239.8/bornesup !nm to eV
bornesup = 1239.8/tempo    !nm to eV


if (type .eq. 'G') then
 write(6,'(A35)') 'Convoluting spectrum using gaussian'
else if (type .eq. 'L') then
 write(6,'(A37)') 'Convoluting spectrum using lorentzian'
endif


write(6,'(a1)') " "

end
