subroutine deal

use declare
implicit none

!Reading input
if (allocated(namegaus)) deallocate(namegaus)
if (allocated(stateresid)) deallocate(stateresid)
if (allocated(Natom)) deallocate(Natom)
if (allocated(res_num)) deallocate(res_num)

!Residus and atoms properties
if (allocated(ground_energy)) deallocate(ground_energy)
if (allocated(energy)) deallocate(energy)
if (allocated(atom_coord)) deallocate(atom_coord)
if (allocated(atom_charge)) deallocate(atom_charge)
if (allocated(masse)) deallocate(masse)
if (allocated(massecenter)) deallocate(massecenter)
if (allocated(coord)) deallocate(coord)
if (allocated(R)) deallocate(R)
if (allocated(Rn)) deallocate(Rn)
if (allocated(el_dip)) deallocate(el_dip)
if (allocated(mag_dip)) deallocate(mag_dip)
if (allocated(rotatory)) deallocate(rotatory)
if (allocated(freq)) deallocate(freq)

!option
if (allocated(comb)) deallocate(comb)

!Hamiltonian
if (allocated(H)) deallocate(H)
if (allocated(eigval)) deallocate(eigval)
if (allocated(work)) deallocate(work)

!Convolution
if (allocated(energy_c)) deallocate(energy_c)
if (allocated(rot_c)) deallocate(rot_c)
if (allocated(lambda_c)) deallocate(lambda_c)
if (allocated(freq_c)) deallocate(freq_c)
if (allocated(absorb)) deallocate(absorb)

!Geometry
if (allocated(at_arom_1)) deallocate(at_arom_1)
if (allocated(at_arom_2)) deallocate(at_arom_2)
if (allocated(normal_vector)) deallocate(normal_vector)
if (allocated(angle)) deallocate(angle)
if (allocated(res_num_ss)) deallocate(res_num_ss)
if (allocated(massecenter_ss)) deallocate(massecenter_ss)
if (allocated(list_atom_ss)) deallocate(list_atom_ss)

if (allocated(list_atom)) deallocate(list_atom)

if (allocated(min_distance_list)) deallocate(min_distance_list)

end
