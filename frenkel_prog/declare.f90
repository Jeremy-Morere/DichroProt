module declare

implicit none

!General
integer :: a,i,j,k,l,m,n,p,s,at
integer :: n_start, n_end
integer :: error,info,lwork
character*64 :: arg1,arg2
character*64 :: output
character*128 :: line

!Reading file
character*64,allocatable :: namegaus(:)
character*64 :: pdb_file, prmtop_file
integer,allocatable :: stateresid(:),Natom(:),res_num(:)
real*8 :: tempo

!Amber/Gaussian
character*64 :: folder_name,res1,res2


!States, residus, atoms
integer :: state,Maxstate,Limstate
integer :: Nresid,Narom,Nss
integer :: atom,Maxatom,AtomNum 

!Residus and atoms properties
real*8,allocatable :: energy(:,:),ground_energy(:)
real*8,allocatable :: el_dip(:,:,:),mag_dip(:,:,:),rotatory(:,:),freq(:,:)
real*8 :: Mass
character*2 :: at_type

!option
logical :: stronghest,extend,allcomb,nocoupling
integer,allocatable :: comb(:)

!Hamiltonian
real*8,allocatable :: H(:,:),eigval(:),work(:)
integer :: d,od1,od2
integer :: r1,r2,s1,s2

!Rotatory and oscillator strenght
real*8 :: PS1,PS2,PS3
real*8 :: cross(3), vect(3)

!Convolution
integer :: npoint,npeak
character*128 :: type
real*8 :: dw,c,q
real*8 :: fwhm,bornesup,borneinf
real*8,allocatable :: energy_c(:),rot_c(:),freq_c(:)
real*8,allocatable :: lambda_c(:),absorb(:)
real*8 :: max_abs

!Geometry
real*8,allocatable :: atom_coord(:,:,:), atom_mass(:,:)
real*8,allocatable :: masse(:), massecenter(:,:)
real*8,allocatable :: coord(:,:)
real*8,allocatable :: R(:,:,:), Rn(:,:)
real*8,allocatable :: at_arom_1(:,:), at_arom_2(:,:)
character*3 :: type_tempo, type_res
character*128 :: list_at,at_num
character*4 :: at_name
integer,allocatable :: res_arom_num(:), res_ss_num(:)
integer :: Nat, Maxat

real*8,allocatable :: normal_vector(:,:)
real*8 :: vect_1(3), vect_2(3)
real*8 :: norm

real*8,allocatable :: angle(:,:,:)

character*256,allocatable :: list_atom(:)

real*8,allocatable:: min_distance_list(:,:)
real*8 :: min_val

!Parameters
real*8,parameter :: pi = 2*dacos(0.0d0)
real*8,parameter :: Na = 6.02214076d23 !mol-1
real*8,parameter :: beta = 1.2563d0

end module declare
