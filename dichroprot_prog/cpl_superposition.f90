subroutine cpl_superposition

!Superposition of CD spectrum without any coupling

use declare
implicit none

write(6,'(A24)') "#-Spectrum superposition"
output="superposition"


!Allocates memory for energy and rotatory strength during convolution of spectra.
allocate(energy_c(Nband))
allocate(energy_tempo(Nband))
allocate(rot_c(Nband))

do ti=1,Nband
 !Residues and states involved in the transition ti.
 ri=transi_band(ti,1)
 si=transi_band(ti,2)

 energy_c(ti) = energy(ri,si)*2.7211407953d1 !au to eV
 rot_c(ti)    = rotatory(ri,si) 
enddo

energy_tempo = energy_c

!Print rotatory strength and each contribution
write(6,'(a2,3a8,2a12,a8,a12,a8,a12,a8,a12,a8)') '| ','residue','state','Energy','Rot'

do i=1,Nband
 ti = MINLOC(energy_tempo,dim=1,mask= energy_tempo > 0)

 energy_tempo(ti) = -1.0d0
 ri=transi_band(ti,1)
 si=transi_band(ti,2)

 if (ri .le. Narom) then
  write(6,'(a2,2i8,f8.2,e12.3)') &
          '| ', res_num(ri), si, 1.239d3/energy_c(ti),rot_c(ti) !Converte into nm
 else
  write(6,'(a2,2i4,i8,f8.2,e12.3)') &
          '| ', res_num(ri), res_num(ri+1), si, 1.239d3/energy_c(ti), rot_c(ti)
 endif
enddo


call convolution

deallocate(energy_c,rot_c,energy_tempo)

end

