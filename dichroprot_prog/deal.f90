subroutine deal

use declare
implicit none

!Reading input
if (allocated(namegaus)) deallocate(namegaus)
if (allocated(stateresid)) deallocate(stateresid)
if (allocated(Natom)) deallocate(Natom)
if (allocated(res_num)) deallocate(res_num)
if (allocated(transi_band)) deallocate(transi_band)


!Residue and atom properties
if (allocated(ground_energy)) deallocate(ground_energy)
if (allocated(energy)) deallocate(energy)
if (allocated(atom_coord)) deallocate(atom_coord)
if (allocated(atom_charge)) deallocate(atom_charge)
if (allocated(chargecenter)) deallocate(chargecenter)
if (allocated(R)) deallocate(R)
if (allocated(Rn)) deallocate(Rn)
if (allocated(el_dip)) deallocate(el_dip)
if (allocated(el_dip_perm)) deallocate(el_dip_perm)
if (allocated(mag_dip)) deallocate(mag_dip)
if (allocated(rotatory)) deallocate(rotatory)
if (allocated(freq)) deallocate(freq)

!Hamiltonian
if (allocated(H)) deallocate(H)
if (allocated(eigval)) deallocate(eigval)
if (allocated(work)) deallocate(work)

!Rotatory and oscillator strengths
if (allocated(Tino)) deallocate(Tino)
if (allocated(Tino_a)) deallocate(Tino_a)
if (allocated(Tino_e)) deallocate(Tino_e)
if (allocated(Tino_b)) deallocate(Tino_b)
if (allocated(Tino_f)) deallocate(Tino_f)
if (allocated(energy_tempo)) deallocate(energy_tempo)
if (allocated(rot_mat)) deallocate(rot_mat)
if (allocated(rot_pert)) deallocate(rot_pert)

!Convolution
if (allocated(energy_c)) deallocate(energy_c)
if (allocated(rot_c)) deallocate(rot_c)
if (allocated(lambda_c)) deallocate(lambda_c)
if (allocated(freq_c)) deallocate(freq_c)
if (allocated(absorb)) deallocate(absorb)
if (allocated(absorb_UV)) deallocate(absorb_UV)

!Geometry
if (allocated(at_arom_1)) deallocate(at_arom_1)
if (allocated(at_arom_2)) deallocate(at_arom_2)
if (allocated(res_arom_num)) deallocate(res_arom_num)
if (allocated(res_ss_num)) deallocate(res_ss_num)
if (allocated(normal_vector)) deallocate(normal_vector)
if (allocated(angle)) deallocate(angle)

if (allocated(list_atom)) deallocate(list_atom)

if (allocated(min_distance_list)) deallocate(min_distance_list)

end
