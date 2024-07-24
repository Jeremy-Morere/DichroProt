subroutine rlog

!Read gaussian 09 .log file.

use declare
implicit none

!FIRST PART: Read numbers of atom and center of charge for each residue
allocate(Natom(Nresid))
Natom = 0

do i=1,Nresid
open(10,file=TRIM(namegaus(i))//"-S1.log",form='formatted')
 do
  read(10,'(A128)',iostat=error) line
  if (error .ne. 0) exit
  if (line(1:8) .eq. ' NAtoms=') then
   read(line,'(10x,i3)') Natom(i)
   exit
  endif
 enddo
close(10)
enddo


!SECONDE PART: Acquiring the data from the gaussian output files.
Maxatom = maxval(Natom)

allocate(energy(Nresid,Nstate))
allocate(ground_energy(Nresid))
allocate(el_dip(Nresid,Nstate,3))
allocate(mag_dip(Nresid,Nstate,3))
allocate(rotatory(Nresid,Nstate))
allocate(freq(Nresid,Nstate))

energy = 0.0d0
el_dip = 0.0d0
el_dip = 0.0d0
rotatory = 0.0d0
freq = 0.0d0
Nband = 0


do ri=1,Nresid
   filename="-S1.log"
   open(10,file=trim(namegaus(ri))//"-S1.log",form='formatted')
   n = 0
   do
      read(10,'(A128)',iostat=error) line
      if (error .ne. 0) exit

      !Read the ground state energy
      if(line(1:10) .eq. ' SCF Done:') then
         read(line,'(23x,f16.9)') ground_energy(ri) !au
      endif

      !Read the electric transition dipole moments
      if (line(1:65) .eq. ' Ground to excited state transition electric dipole moments (Au):' ) then
         read(10,'(A128)') line
         do si=1,Nstate
            read(10,'(13x,3f12.4)') el_dip(ri,si,1),el_dip(ri,si,2),el_dip(ri,si,3)
         enddo
      endif

      !Read magnetic transition dipole moments
      if (line(1:65) .eq. ' Ground to excited state transition magnetic dipole moments (Au):' ) then
         read(10,'(A128)') line
         do si=1,Nstate
            read(10,'(13x,3f12.4)') mag_dip(ri,si,1),mag_dip(ri,si,2),mag_dip(ri,si,3)
         enddo
      endif

      !Read the transition energies and oscillator strengths
      if (line(1:14) .eq. ' Excited State')  then
         n = n + 1
         read(line,'(40x,f7.4,18x,f6.4)') energy(ri,n), freq(ri,n)
         energy(ri,n) = energy(ri,n)/2.72116d1 !eV to au
         !Count if transition in near UV
         if (energy(ri,n) .lt. 1.981020d-1) then
          Nband = Nband + 1
         endif
      endif

      !Read  the oscillator strenghts
      if (line(47:62) .eq. 'ZZ     R(length)') then
         do si=1,Nstate
            read(10,'(A128)') line
            read(line,'(52x,f9.4)') rotatory(ri,si) !au
         enddo
      endif
   enddo
enddo

close(10)


!Read permanant electric dipole 
if (full) then
 allocate(el_dip_perm(Nresid,0:Nstate,3))
 el_dip_perm = 0.0d0
 
 do ri=1,Nresid
  do si=0,Nstate
 
   !Generate gaussian log name
   if (si .lt. 10) then
    write(filename,'(a2,i1,a4)') "-S", si, ".log"
   else
    write(filename,'(a2,i2,a4)') "-S", si, ".log"
   endif
 
   filename = trim(namegaus(ri))//filename
   
   !Open file and read permanent dipole.
   open(10,file=filename,form='formatted')
    do
     read(10,'(A128)',iostat=error) line
     if (error .ne. 0) exit
 
     if (line(1:14) .eq. " Dipole moment") then
      read(10,'(A128)',iostat=error) line
      read(line,'(6x,f20.4,6x,f20.4,6x,f20.4)') el_dip_perm(ri,si,1), el_dip_perm(ri,si,2), el_dip_perm(ri,si,3) !in debye
      exit
     endif
   enddo
 
   close(10)
  enddo
 enddo
 
 el_dip_perm = el_dip_perm*3.93456d-1 !debye to au
endif

!Select transition in spectral band
allocate(transi_band(Nband,2))

write(6,'(a1)') " "
write(6,'(a22)') "Transitions in near-UV"
write(6,'(a18)') "Residu       State"

n = 0
do ri=1,Nresid
 do si=1,Nstate
  if (energy(ri,si) .lt. 1.981020d-1) then
   n = n + 1
   transi_band(n,1) = ri
   transi_band(n,2) = si

   if (ri .le. Narom) then
    write(6,'(i6,6x,i6)') res_num(ri), si
   else
    write(6,'(3i6)') res_num(ri), res_num(ri+1), si
   endif
  endif
 enddo
enddo


!zone test

if (.false.) then
 write(6,'(a1)') " "
 open(10,file="UV.dat",form='formatted')
 do ti=1,Nband
  ri=transi_band(ti,1)
  si=transi_band(ti,2)
  write(10,'(f8.4,f8.4)') 45.5634/energy(ri,si), freq(ri,si)
 enddo
 close(10)
endif




if (.false.) then
do ti=1,Nband
 ri=transi_band(ti,1)
 si=transi_band(ti,2)
write(*,*) "############################################"
write(*,'(a64)') namegaus(ri)
write(6,'(i4,a3,2i4)') ti," : ",ri,si
  
write(*,'(a6)') "el dip"
  write(*,'(3f12.4)') el_dip(ri,si,1),el_dip(ri,si,2),el_dip(ri,si,3)
write(*,'(a7)')   "mag dip"
  write(*,'(3f12.4)') mag_dip(ri,si,1),mag_dip(ri,si,2),mag_dip(ri,si,3)
write(*,'(a20)')    "energy de transition"
  write(*,'(f12.4)') energy(ri,si)*2.72116d1 
if (full) then
write(*,'(a15)')    "dipole pemanent"
  write(*,'(3f12.4)') el_dip_perm(ri,si,1)/3.93456d-1, el_dip_perm(ri,si,2)/3.93456d-1, el_dip_perm(ri,si,3)/3.93456d-1
endif
write(*,*)
enddo
endif

deallocate(namegaus)

end
