subroutine rlog

!Read gaussian 09 .log file.

use declare
implicit none

!FIRST PART: Find numbers of transition and atom for each residue
allocate(stateresid(Nresid))
allocate(Natom(Nresid))

do i=1,Nresid
   n = 0
   open(10,file=namegaus(i),form='formatted')
   do
      read(10,'(A128)',iostat=error) line
      if (error .ne. 0) exit

      if (line(1:8) .eq. ' NAtoms=') then
         read(line,'(10x,i3)') Natom(i)

         !Substract H atom add for qmmm
         if (i .le. Narom) Natom(i) = Natom(i) - 1
         if (i .gt. Narom) Natom(i) = Natom(i) - 2
      endif

      if (line(37:65) .eq.  'electric dipole moments (Au):') then
         do
            read(10,'(A128)') line
            if (line(37:65) .eq. 'velocity dipole moments (Au):') exit
            n = n + 1
         enddo
         stateresid(i) = n - 1
      endif
   enddo
close(10)
enddo


Maxstate = maxval(stateresid)
Maxatom = maxval(Natom)

allocate(energy(Nresid,Maxstate))
allocate(ground_energy(Nresid))
allocate(el_dip(Nresid,Maxstate,3))
allocate(rotatory(Nresid,Maxstate))
allocate(freq(Nresid,Maxstate))

energy = 0.0d0
el_dip = 0.0d0
rotatory = 0.0d0
freq = 0.0d0
        
!SECONDE PART: Acquiring the data from the gaussian output files.
do i=1,Nresid
   open(10,file=namegaus(i),form='formatted')
   n = 0
   do
      read(10,'(A128)',iostat=error) line
      if (error .ne. 0) exit

      !Read the ground state energy
      if(line(1:10) .eq. ' SCF Done:') then
         read(line,'(23x,f16.9)') ground_energy(i) !au
      endif

      !Read the electric transition dipole moments
      if (line(1:65) .eq. ' Ground to excited state transition electric dipole moments (Au):' ) then
         read(10,'(A128)') line
         do k=1,stateresid(i)
            read(10,'(13x,3f12.4)') el_dip(i,k,1),el_dip(i,k,2),el_dip(i,k,3)
         enddo
      endif

      !Read the transition energies and oscillator strengths
      if (line(1:14) .eq. ' Excited State')  then
         n = n + 1
         read(line,'(40x,f7.4,18x,f6.4)') energy(i,n), freq(i,n)
         energy(i,n) = energy(i,n)/2.72116d1 !eV to au
      endif

      !Read  the oscillator strenghts
      if (line(47:62) .eq. 'ZZ     R(length)') then
         do k=1,stateresid(i)
            read(10,'(A128)') line
            read(line,'(52x,f9.4)') rotatory(i,k) !au
         enddo
      endif
   enddo
enddo

close(10)

deallocate(namegaus)

end
