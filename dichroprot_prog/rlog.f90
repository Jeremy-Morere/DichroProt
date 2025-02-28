subroutine rlog

!Read Gaussian09 .log file.

use declare
implicit none

!FIRST PART: Find numbers of transition and atom for each residue
allocate(Natom(Nresid))

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

Maxatom = maxval(Natom)

allocate(energy(Nresid,Nstate))
allocate(ground_energy(Nresid))
allocate(el_dip(Nresid,Nstate,3))
allocate(mag_dip(Nresid,Nstate,3))
allocate(rotatory(Nresid,Nstate))
allocate(freq(Nresid,Nstate))

energy   = 0.0d0
el_dip   = 0.0d0
el_dip   = 0.0d0
mag_dip  = 0.0d0
rotatory = 0.0d0
freq     = 0.0d0
Nband = 0        

!SECOND PART: Acquiring the data from the Gaussian09 output files.
do ri=1,Nresid
   filename="-S1.log"
   open(10,file=trim(namegaus(ri))//"-S1.log",form='formatted')
   n = 0
   do
      read(10,'(A128)',iostat=error) line
      if (error .ne. 0) exit

      !Read the electric transition dipole moments
      if (line(1:65) .eq. ' Ground to excited state transition electric dipole moments (Au):' ) then
         read(10,'(A128)') line
         do si=1,Nstate
            read(10,'(13x,3f12.4)') el_dip(ri,si,1),el_dip(ri,si,2),el_dip(ri,si,3) !au
         enddo
      endif

      !Read magnetic transition dipole moments
      if (line(1:65) .eq. ' Ground to excited state transition magnetic dipole moments (Au):' ) then
         read(10,'(A128)') line
         do si=1,Nstate
            read(10,'(13x,3f12.4)') mag_dip(ri,si,1),mag_dip(ri,si,2),mag_dip(ri,si,3) !electron magnetic moment unit
         enddo
      endif

      !Read the transition energies and oscillator strengths
      if (line(1:14) .eq. ' Excited State')  then
         n = n + 1
         read(line,'(40x,f7.4,18x,f6.4)') energy(ri,n), freq(ri,n) !eV, unitless
         energy(ri,n) = energy(ri,n)/2.72116d1 !eV to au
         !Count if transition in near UV
         if (energy(ri,n) .lt. 1.981020d-1) then !<230nm
          Nband = Nband + 1
         endif
      endif

      !Read the oscillator strenghts
      if (line(47:62) .eq. 'ZZ     R(length)') then
         do si=1,Nstate
            read(10,'(A128)') line
            read(line,'(52x,f9.4)') rotatory(ri,si) !cgs
         enddo
      endif
   enddo
enddo
mag_dip = mag_dip/2.00231930436152d0 !electron magnetic momentunit to au
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
 
     if (line(1:48) .eq. " Dipole moment (field-independent basis, Debye):") then
      read(10,'(A128)',iostat=error) line
      read(line,'(6x,f20.4,6x,f20.4,6x,f20.4)') el_dip_perm(ri,si,1), el_dip_perm(ri,si,2), el_dip_perm(ri,si,3) !in Debye
      exit
     endif
   enddo
 
   close(10)
  enddo
 enddo
 el_dip_perm = el_dip_perm*3.93456d-1 !Converts electric dipole moment from Debye to atomic units.
endif

!Select transition in spectral band
allocate(transi_band(Nband,2))
n = 0
do ri=1,Nresid
 do si=1,Nstate
  if (energy(ri,si) .lt. 1.981020d-1) then
   n = n + 1
   transi_band(n,1) = ri
   transi_band(n,2) = si
  endif
 enddo
enddo

deallocate(namegaus)

end
