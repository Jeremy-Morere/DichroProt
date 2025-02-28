The source repository  contains all programs and executables necessary to obtain a near-UV CD spectrum.

dichroprot_prog	: Coupled gaussian 09 calculations; Generates a CD spectrum with or without coupling.
geometry_prog	: Analyses the protein structure from a .pdb; generates a geometry.dat file used to obtain the CD spectrum.
other_prog	: Other small programs for post-processing.
exectutable	: Contains all executables necessary to scan a trajectory (see Tutorial for more detail)
test_snapshot	: Test repository (Lysozyme (PDB:2YVB) is given in example)

#--Prerequisite--#
The current version of DichroProt has been developed ans tested under Centos Linux 7 with 
	mkl libraries: 2019.4
	Fortran compiler : Intel 2019.4
	Amber : 16
	Gaussian : 09

Paths to mkl libraries, Fortran compiler, Gaussian and Amber has to be precise in the install file.

The executable file has to be add to PATH. It can be done by adding to the .bashrc file
	export PATH=/path_to_executable:$PATH


#--Tutorial--#

Necessary input files for a scan:
        - topology file (.prmtop)
        - trajectory (.dcd or .nc)
Necessary input files for a snapshot
        - topology file (.prmtop)
        - geometry file (.pdb)

--Description of executables--
geometry:	 Spot Photo-active residus (PAR) and compute center of charge for a given conformation and other information.
			-Inputs: .pdb file
			-Outputs: geometry.dat file
create_input:	 Generate global input, can be use for scan, qmmm or calculation.
                        -Outputs: GlobalInput
			-call: create_input

prep_scan:	 Prepare repositories for each conformation of a trajectory.
                        -Inputs: .prmtop and .dcd/.nc files
                        -Outputs: one repository per conformation (step of 20 along the trajectory) with every useful file for prep_qmmm exectutable
                        -call: prep_scan GlobalInput
prep_qmmm:	 Prepare repositories for each state for a given conformation (need: .prmtop and .pdb files ; call: prep_qmmm QmmmInput)
                        -Inputs: .prmtop and .pdb files
                        -Outputs: one repository per PAR, one calculation repository per requested state and QmmmInput
                        -call: prep_qmmm QmmmInput
prep_dichroprot: Gather G09 output in a Dichroprot repository.
                        -need: launch_qmmm output
                        -Outputs: Dichroprot calculation repository and DichroProtInput
                        -call: prep_dichroprot QmmmInput

launch_scan:	 launch QMMM calculation of a trajectory.
                        -Inputs: prep_scan output
                        -call: launch_scan GlobalInput
launch_qmmm:	 launch QMMM calculation for a given conformation.
                        -Inputs: prep_qmmm output
                        -call: launch_qmmm QmmmInput
launch_dichroprot: launch CD calculation of a trajectory.
                        -Inputs: launch_scan output
                        -call: launch_dichroprot GlobalInput

dichroprot:	 Compute stick spectrum and convoluted spectrum for a given conformation.
                        -call: dichroprot DichroProtInput


#--Scan tutorial--#
1) Enter your working repository containing trajectory and topology files (.prmtop and .dcd/.nc)
2) Launch `create_input` and adapt the GlobalInput file to your portein
        -PRMTOP:       		protein.prmtop
        -trajectory:   		protein.nc
        -Level of theory: 	B3LYP/6-311+G**
        -Number of state:	5	(Advice: A lot of transitions are missed below 5)
        -Full contribution: 	n	(y : protein with few PAR or highly interacting PAR)
3) Launch `prep_scan GlobalInput`
4) Adapt the submith.sh to your architecture and gau_job.tpl if you need to add specific keywords for the gaussian calculation.
5) Launch `launch_scan GlobalInput` (It will launch every qmmm calculations)
6) When every Qmmm calculations end, check if any crashed and relaunch them, following step 1-2 of Conformation tutorial.
7) Launch `launch_dichroprot GlobalInput`
8) iAll output files will be gathered in the spectrum repository.


#--Snapshot tutorial--#
	The lysozymz protein (PDB:2YVB) is given in example in the snapshot_test repository 
1) Enter your working repository containing structure and topology files (.prmtop and .pdb)
2) Launch `create_input` and adapt GlobalInput file to your portein
	-PDB :           protein.pdb
        -PRMTOP :        protein.prmtop
        -Level of theory : B3LYP/6-311+G**
        -Number of state : 5                    (Advice: A lot of transitions are missed below 5)
        -Full contribution : n                  (y : protein with few PAR or highly interacting PAR)
3) Launch `prep_qmmm GlobalInput`
4) Adapt the submith.sh to your architecture and gau_job.tpl if you need to add specific keywords for the gaussian calculation.
5) Launch `launch_qmmm GlobalInput`
7) Launch `prep_dichroprot GlobalInput`
8) All output files will be gathered in the DichroProt repository.
9) Launch 'dichroprot DichroProtInput'


