module declare

implicit none

!General
integer       :: a,i,j,k,l,m,n,p,s,at
integer       :: n_start, n_end
integer       :: error,info,lwork
character*64  :: arg1,arg2
character*64  :: output
character*128 :: line

!Reading file
character*64,allocatable :: namegaus(:)
character*64             :: pdb_file,filename
integer,allocatable      :: stateresid(:),Natom(:),res_num(:)
real*8                   :: tempo

!States, residus, atoms
integer             :: state,Maxstate,Nstate
integer             :: Nresid,Narom,Nss
integer             :: atom,Maxatom,AtomNum 
integer             :: Nband
integer,allocatable :: transi_band(:,:)

!Residue and atom properties
real*8,allocatable :: energy(:,:), ground_energy(:), rotatory(:,:), freq(:,:)
real*8,allocatable :: el_dip_perm(:,:,:), el_dip(:,:,:), mag_dip(:,:,:)
real*8             :: Charge
character*2        :: at_type

!option
character*1         :: matrix_t, perturbation_t, superposition_t, full_t
logical             :: matrix, perturbation, superposition, full

!Hamiltonian
real*8,allocatable :: H(:,:),eigval(:),work(:)
integer            :: d,od1,od2
integer            :: r1,r2,ri,s1,s2,si,t1,t2,ti

!Rotatory and oscillator strengths
real*8             :: cross(3), vect(3)
real*8,allocatable :: Tino(:), Tino_a(:), Tino_e(:), Tino_b(:), Tino_f(:)
real*8             :: mat, mat_mag_el, mat_el_el
real*8,allocatable :: energy_tempo(:)
real*8,allocatable :: rot_mat(:), rot_pert(:)

!Convolution
integer            :: npoint,npeak
character*128      :: type
real*8             :: dw,c,q
real*8             :: fwhm, bornesup, borneinf
real*8,allocatable :: energy_c(:), rot_c(:), freq_c(:)
real*8,allocatable :: lambda_c(:), absorb(:), absorb_UV(:)
real*8             :: max_abs, min_abs, coeff

!Geometry
real*8,allocatable  :: atom_coord(:,:,:), atom_charge(:,:)
real*8,allocatable  :: chargecenter(:,:)
real*8,allocatable  :: R(:,:,:), Rn(:,:)
real*8,allocatable  :: at_arom_1(:,:), at_arom_2(:,:)
character*3         :: type_tempo, type_res
character*128       :: list_at,at_num
character*4         :: at_name
integer,allocatable :: res_arom_num(:), res_ss_num(:)
integer             :: Nat, Maxat

real*8,allocatable :: normal_vector(:,:)
real*8             :: vect_1(3), vect_2(3)
real*8             :: norm

real*8,allocatable :: angle(:,:,:)

character*256,allocatable :: list_atom(:)

real*8,allocatable:: min_distance_list(:,:)
real*8            :: min_val

!Parameters
real*8,parameter :: pi = dacos(-1.0d0)
real*8,parameter :: Na = 6.02214076d23 !mol-1
real*8,parameter :: beta = 1.2563d0
real*8,parameter :: au_cgs = 4.714436d2

!--forbiden name--!
! res1, s1i, s1f, res2, s2i, s2f
! PS1, PS2, PS3
! rot

end module declare
