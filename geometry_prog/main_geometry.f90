program main_geometry

!Get the coordinate of each photo-active side chain

use declare
implicit none
real*8 :: cnorm

call getarg(1,pdb_file)

!STEP 1: Find number of Photo-Active Residue
call find_number_PAR

!STEP 2: Get atoms coordinateis and number
call rpdb

!STEP 3: Compute center of molecule
call coc

!STEP 4: find partner of each cystein
call store_ss

!STEP 5: Select residu to put in the same qmmm partition

call woutput


call deal

end
